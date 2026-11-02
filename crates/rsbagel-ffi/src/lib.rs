//! Thin in-memory C ABI. Ownership and buffer layouts: docs/screenforge_integration.md.
//! BAGEL2-derived behavior: Copyright (c) 2020 Hart Lab, MIT; see docs/BAGEL2-LICENSE.txt.

use rsbagel_core::{BfConfig, Error, FcConfig, Matrix};
use serde_json::{json, Value};
use std::ffi::c_void;
use std::mem::{align_of, size_of};
use std::ops::ControlFlow;
use std::panic::{catch_unwind, AssertUnwindSafe};
use std::ptr;
use std::slice;
use std::sync::atomic::{AtomicBool, Ordering};

pub const RSB_OK: i32 = 0;
pub const RSB_INVALID_ARGUMENT: i32 = 1;
pub const RSB_NUMERICAL_ERROR: i32 = 2;
pub const RSB_LEGACY_FAILURE: i32 = 3;
pub const RSB_CANCELLED: i32 = 4;
pub const RSB_IO_ERROR: i32 = 5;
pub const RSB_INTERNAL_ERROR: i32 = 6;

#[repr(C)]
#[derive(Clone, Copy)]
pub struct RsBytes {
    pub data: *const u8,
    pub len: usize,
}
#[repr(C)]
#[derive(Clone, Copy)]
pub struct RsStrings {
    pub data: *const RsBytes,
    pub len: usize,
}
#[repr(C)]
#[derive(Clone, Copy)]
pub struct RsF64 {
    pub data: *const f64,
    pub len: usize,
}
/// Labels determine rows/columns; `values` is a flat row-major buffer.
#[repr(C)]
#[derive(Clone, Copy)]
pub struct RsMatrix {
    pub row_ids: RsStrings,
    pub row_genes: RsStrings,
    pub columns: RsStrings,
    pub values: RsF64,
}
#[repr(C)]
#[derive(Clone, Copy)]
pub struct RsNumericView {
    pub data: *const f64,
    pub len: usize,
    pub rows: usize,
    pub columns: usize,
}

struct NumericBuffer {
    name: &'static str,
    rows: usize,
    columns: usize,
    values: Vec<f64>,
}
/// Opaque immutable result, including errors. Release with `rsbagel_result_free`.
pub struct RsResult {
    status: i32,
    error: Vec<u8>,
    metadata: Vec<u8>,
    buffers: Vec<NumericBuffer>,
}
/// A cancellation token can be requested concurrently, but must outlive its run.
pub struct RsCancelToken {
    cancelled: AtomicBool,
}
pub type RsProgress = Option<unsafe extern "C" fn(*mut c_void, usize, usize) -> i32>;

type FfiResult<T> = Result<T, FfiError>;
struct FfiError {
    status: i32,
    message: String,
}
impl FfiError {
    fn invalid(message: impl Into<String>) -> Self {
        Self {
            status: RSB_INVALID_ARGUMENT,
            message: message.into(),
        }
    }
    fn internal(message: impl Into<String>) -> Self {
        Self {
            status: RSB_INTERNAL_ERROR,
            message: message.into(),
        }
    }
}
impl From<Error> for FfiError {
    fn from(error: Error) -> Self {
        let status = match &error {
            Error::Invalid(_) | Error::Csv(_) | Error::Json(_) => RSB_INVALID_ARGUMENT,
            Error::Numerical(_) => RSB_NUMERICAL_ERROR,
            Error::Legacy(_) => RSB_LEGACY_FAILURE,
            Error::Cancelled => RSB_CANCELLED,
            Error::Io(_) => RSB_IO_ERROR,
        };
        Self {
            status,
            message: error.to_string(),
        }
    }
}

/// Checks representation constraints before creating a borrowed slice. The C
/// caller remains responsible for the actual allocation and its live lifetime.
fn check_pointer<T>(data: *const T, len: usize, label: &str) -> FfiResult<()> {
    if len == 0 {
        return Ok(());
    }
    if data.is_null() {
        return Err(FfiError::invalid(format!(
            "{label}: null pointer with nonzero length"
        )));
    }
    if !(data as usize).is_multiple_of(align_of::<T>()) {
        return Err(FfiError::invalid(format!(
            "{label}: pointer alignment is invalid"
        )));
    }
    let bytes = len
        .checked_mul(size_of::<T>())
        .ok_or_else(|| FfiError::invalid(format!("{label}: byte length overflow")))?;
    if bytes > isize::MAX as usize || (data as usize).checked_add(bytes).is_none() {
        return Err(FfiError::invalid(format!(
            "{label}: slice exceeds addressable range"
        )));
    }
    Ok(())
}
unsafe fn checked_slice<'a, T>(data: *const T, len: usize, label: &str) -> FfiResult<&'a [T]> {
    check_pointer(data, len, label)?;
    if len == 0 {
        Ok(&[])
    } else {
        Ok(slice::from_raw_parts(data, len))
    }
}
unsafe fn strings(view: RsStrings, label: &str) -> FfiResult<Vec<String>> {
    let views = checked_slice(view.data, view.len, label)?;
    views
        .iter()
        .enumerate()
        .map(|(index, value)| {
            let bytes = checked_slice(value.data, value.len, label)?;
            std::str::from_utf8(bytes)
                .map(str::to_owned)
                .map_err(|_| FfiError::invalid(format!("{label}[{index}]: invalid UTF-8")))
        })
        .collect()
}
unsafe fn matrix(view: RsMatrix) -> FfiResult<Matrix> {
    let _profile_span = rsbagel_core::profiling::span("ffi_matrix_conversion");
    let rows = view.row_ids.len;
    let columns = view.columns.len;
    let count = rows
        .checked_mul(columns)
        .ok_or_else(|| FfiError::invalid("matrix dimensions overflow"))?;
    if view.row_genes.len != rows || view.values.len != count {
        return Err(FfiError::invalid(
            "matrix dimensions do not match buffer lengths",
        ));
    }
    // Check every outer view before reading the labels or constructing slices.
    check_pointer(view.row_ids.data, rows, "row_ids")?;
    check_pointer(view.row_genes.data, rows, "row_genes")?;
    check_pointer(view.columns.data, columns, "columns")?;
    check_pointer(view.values.data, count, "values")?;
    let result = Matrix {
        index_name: "REAGENT_ID".into(),
        row_ids: strings(view.row_ids, "row_ids")?,
        row_genes: strings(view.row_genes, "row_genes")?,
        columns: strings(view.columns, "columns")?,
        values: checked_slice(view.values.data, count, "values")?.to_vec(),
    };
    result.validate()?;
    Ok(result)
}
unsafe fn json_config(view: RsBytes) -> FfiResult<Value> {
    let bytes = checked_slice(view.data, view.len, "config JSON")?;
    if bytes.is_empty() {
        return Ok(json!({}));
    }
    serde_json::from_slice(bytes)
        .map_err(|e| FfiError::invalid(format!("invalid config JSON: {e}")))
}
fn fc_config(value: Value, input: &Matrix) -> FfiResult<FcConfig> {
    let object = value
        .as_object()
        .ok_or_else(|| FfiError::invalid("FC config must be a JSON object"))?;
    for key in object.keys() {
        if !matches!(
            key.as_str(),
            "controls" | "control_labels" | "min_reads" | "pseudocount"
        ) {
            return Err(FfiError::invalid(format!("unknown FC config field: {key}")));
        }
    }
    let controls: Vec<usize> = match object.get("controls") {
        Some(value) => serde_json::from_value(value.clone())
            .map_err(|e| FfiError::invalid(format!("invalid controls: {e}")))?,
        None => vec![0],
    };
    let control_labels: Vec<String> = match object.get("control_labels") {
        Some(value) => serde_json::from_value(value.clone())
            .map_err(|e| FfiError::invalid(format!("invalid control_labels: {e}")))?,
        None => controls
            .iter()
            .map(|&i| {
                input
                    .columns
                    .get(i)
                    .cloned()
                    .ok_or_else(|| FfiError::invalid("FC control column out of range"))
            })
            .collect::<FfiResult<_>>()?,
    };
    let number = |key: &str, default: f64| -> FfiResult<f64> {
        match object.get(key) {
            None => Ok(default),
            Some(value) => value
                .as_f64()
                .ok_or_else(|| FfiError::invalid(format!("{key} must be a number"))),
        }
    };
    Ok(FcConfig {
        controls,
        control_labels,
        min_reads: number("min_reads", 0.0)?,
        pseudocount: number("pseudocount", 5.0)?,
    })
}
fn buffer(
    name: &'static str,
    rows: usize,
    columns: usize,
    values: Vec<f64>,
) -> FfiResult<NumericBuffer> {
    if rows.checked_mul(columns) != Some(values.len()) {
        return Err(FfiError::internal(format!(
            "core result {name} has inconsistent dimensions"
        )));
    }
    Ok(NumericBuffer {
        name,
        rows,
        columns,
        values,
    })
}
fn matrix_metadata(matrix: &Matrix) -> Value {
    json!({"index_name":matrix.index_name,"row_ids":matrix.row_ids,"row_genes":matrix.row_genes,"columns":matrix.columns})
}
fn success(mut metadata: Value, buffers: Vec<NumericBuffer>) -> FfiResult<RsResult> {
    metadata["abi_version"] = json!(1);
    metadata["source_commit"] = json!(rsbagel_core::SOURCE_COMMIT);
    metadata["source_dirty"] = json!(rsbagel_core::SOURCE_DIRTY);
    metadata["ffi_version"] = json!(env!("CARGO_PKG_VERSION"));
    metadata["numeric_order"] = json!("row-major f64");
    metadata["nonfinite_policy"] = json!("Numeric buffers preserve IEEE-754 values. JSON nonfinite values are null; optional-field presence is recorded separately.");
    metadata["buffers"] = Value::Array(buffers.iter().map(|b| json!({"name":b.name,"rows":b.rows,"columns":b.columns,"length":b.values.len()})).collect());
    let metadata = serde_json::to_vec(&metadata).map_err(|e| FfiError::internal(e.to_string()))?;
    Ok(RsResult {
        status: RSB_OK,
        error: vec![],
        metadata,
        buffers,
    })
}
fn failure(error: FfiError) -> RsResult {
    let metadata = json!({"abi_version":1,"status":error.status,"error":error.message});
    RsResult {
        status: error.status,
        error: error.message.into_bytes(),
        metadata: metadata.to_string().into_bytes(),
        buffers: vec![],
    }
}
fn execute(action: impl FnOnce() -> FfiResult<RsResult>) -> *mut RsResult {
    let result = match catch_unwind(AssertUnwindSafe(action)) {
        Ok(Ok(result)) => result,
        Ok(Err(error)) => failure(error),
        Err(payload) => {
            let message = payload
                .downcast_ref::<String>()
                .cloned()
                .or_else(|| payload.downcast_ref::<&str>().map(|s| (*s).to_owned()))
                .unwrap_or_else(|| "Rust panic".into());
            failure(FfiError::internal(format!(
                "contained Rust panic: {message}"
            )))
        }
    };
    Box::into_raw(Box::new(result))
}

#[no_mangle]
pub extern "C" fn rsbagel_abi_version() -> u32 {
    1
}

/// # Safety
/// All nonempty views must refer to live, aligned allocations for the duration
/// of this synchronous call. See the C header for lengths and ownership.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_fc(input: RsMatrix, config_json: RsBytes) -> *mut RsResult {
    let _profile_session = rsbagel_core::profiling::Session::new();
    let _profile_span = rsbagel_core::profiling::span("ffi_workflow");
    execute(|| {
        let input = matrix(input)?;
        let config = fc_config(json_config(config_json)?, &input)?;
        let result = rsbagel_core::run_fc(&input, &config)?;
        let _profile_conversion = rsbagel_core::profiling::span("ffi_fc_result_conversion");
        let cli_input = rsbagel_core::fc_to_bf_input(&result.fold_change);
        let metadata = json!({"workflow":"fc","fold_change":matrix_metadata(&result.fold_change),"normalized":matrix_metadata(&result.normalized),"fold_change_buffer_stage":"raw before TSV rounding","cli_bf_input_available":cli_input.is_ok(),"cli_bf_input_error":cli_input.as_ref().err().map(ToString::to_string),"cli_bf_input_precision":3});
        let mut buffers = vec![
            buffer(
                "fold_change",
                result.fold_change.rows(),
                result.fold_change.cols(),
                result.fold_change.values,
            )?,
            buffer(
                "normalized",
                result.normalized.rows(),
                result.normalized.cols(),
                result.normalized.values,
            )?,
        ];
        if let Ok(input) = cli_input {
            buffers.push(buffer(
                "cli_bf_input",
                input.rows(),
                input.cols(),
                input.values,
            )?);
        }
        success(metadata, buffers)
    })
}

/// # Safety
/// Views/tokens must remain live throughout the call. A progress callback may
/// run on a worker thread and must be thread-safe, must not unwind, and must
/// not free inputs, the token, or an in-flight result. Return nonzero to cancel.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_bf(
    input: RsMatrix,
    essential: RsStrings,
    nonessential: RsStrings,
    config_json: RsBytes,
    callback: RsProgress,
    context: *mut c_void,
    cancel_token: *const RsCancelToken,
) -> *mut RsResult {
    let _profile_session = rsbagel_core::profiling::Session::new();
    let _profile_span = rsbagel_core::profiling::span("ffi_workflow");
    execute(|| {
        let input = matrix(input)?;
        let essential = strings(essential, "essential genes")?;
        let nonessential = strings(nonessential, "nonessential genes")?;
        let config: BfConfig = serde_json::from_value(json_config(config_json)?)
            .map_err(|e| FfiError::invalid(format!("invalid BF config: {e}")))?;
        if !cancel_token.is_null() {
            check_pointer(cancel_token, 1, "cancel token")?;
        }
        let token_address = cancel_token as usize;
        let context_address = context as usize;
        let progress = |completed: usize, total: usize| {
            if token_address != 0
                && (*(token_address as *const RsCancelToken))
                    .cancelled
                    .load(Ordering::Relaxed)
            {
                return ControlFlow::Break(());
            }
            if let Some(function) = callback {
                if function(context_address as *mut c_void, completed, total) != 0 {
                    return ControlFlow::Break(());
                }
            }
            if token_address != 0
                && (*(token_address as *const RsCancelToken))
                    .cancelled
                    .load(Ordering::Relaxed)
            {
                ControlFlow::Break(())
            } else {
                ControlFlow::Continue(())
            }
        };
        if let ControlFlow::Break(()) = progress(0, 0) {
            return Err(Error::Cancelled.into());
        }
        let result =
            rsbagel_core::run_bf(&input, &essential, &nonessential, &config, Some(&progress))?;
        if token_address != 0
            && (*(token_address as *const RsCancelToken))
                .cancelled
                .load(Ordering::Relaxed)
        {
            return Err(Error::Cancelled.into());
        }
        let _profile_conversion = rsbagel_core::profiling::span("ffi_bf_result_conversion");
        let mut gene_values = Vec::with_capacity(result.genes.len() * 7);
        let mut gene_metadata = Vec::with_capacity(result.genes.len());
        for gene in &result.genes {
            gene_values.extend_from_slice(&[
                gene.raw_mean,
                gene.mean,
                gene.std,
                gene.normalized.unwrap_or(f64::NAN),
                gene.boosted.unwrap_or(f64::NAN),
                gene.boosted_std.unwrap_or(f64::NAN),
                gene.boosted_normalized.unwrap_or(f64::NAN),
            ]);
            gene_metadata.push(json!({"gene":gene.gene,"count":gene.count,"normalized_present":gene.normalized.is_some(),"boosted_present":gene.boosted.is_some(),"boosted_std_present":gene.boosted_std.is_some(),"boosted_normalized_present":gene.boosted_normalized.is_some(),"network_evidence":gene.network_evidence,"network_contribution_count":gene.network_contribution_count,"network_zero_applied":gene.network_zero_applied}));
        }
        let replicates = config.columns.len();
        let mut guide_means = Vec::new();
        let mut guide_stds = Vec::new();
        let mut guide_totals = Vec::new();
        let mut guide_metadata = Vec::new();
        for guide in &result.guides {
            if guide.mean_replicates.len() != replicates || guide.std_replicates.len() != replicates
            {
                return Err(FfiError::internal(
                    "core guide replicate dimensions disagree",
                ));
            }
            guide_means.extend_from_slice(&guide.mean_replicates);
            guide_stds.extend_from_slice(&guide.std_replicates);
            guide_totals.push(guide.total);
            guide_metadata.push(json!({"guide":guide.guide,"gene":guide.gene,"count":guide.count}));
        }
        let ranking_column = rsbagel_core::output::default_pr_column(&config);
        let replicate_columns = if config.column_labels.is_empty() {
            config
                .columns
                .iter()
                .map(|&index| input.columns[index].clone())
                .collect::<Vec<_>>()
        } else {
            config.column_labels.clone()
        };
        let replicate_equalization_factor = rsbagel_core::output::replicate_factor(&config);
        let pr_input = rsbagel_core::bf_to_pr_input(&result, &config, ranking_column);
        let metadata = json!({"workflow":"bf","cli_pr_input_available":pr_input.is_ok(),"cli_pr_input_error":pr_input.as_ref().err().map(ToString::to_string),"cli_pr_genes":pr_input.as_ref().ok().map(|v|&v.0),"cli_pr_score_column":ranking_column,"numeric_buffer_stage":"buffers0to3:upstream pre-serialization summaries;cli_pr_input:post-scaling/3decimal/alphabetical","replicate_equalization_factor":replicate_equalization_factor,"serialization_scaling":{"mean":"multiply_factor","boosted":"multiply_factor","boosted_std":"multiply_factor","boosted_normalized":"multiply_factor","std":"unscaled","normalized":"unscaled","guide_totals":"multiply_factor","guide_means":"unscaled","guide_stds":"unscaled"},"genes":gene_metadata,"guides":guide_metadata,"gene_columns":["raw_mean","mean","std","normalized","boosted","boosted_std","boosted_normalized"],"replicate_columns":replicate_columns,"network_profile":config.network_profile,"ranking_column":ranking_column,"diagnostics":result.diagnostics,"correction_coefficients":result.correction_coefficients});
        let mut buffers = vec![
            buffer("genes", result.genes.len(), 7, gene_values)?,
            buffer("guide_means", result.guides.len(), replicates, guide_means)?,
            buffer("guide_stds", result.guides.len(), replicates, guide_stds)?,
            buffer("guide_totals", result.guides.len(), 1, guide_totals)?,
        ];
        if let Ok((genes, scores)) = pr_input {
            buffers.push(buffer("cli_pr_input", genes.len(), 1, scores)?);
        }
        success(metadata, buffers)
    })
}

/// # Safety
/// Each view must satisfy the live-allocation requirements in `rsbagel_fc`.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_pr(
    genes: RsStrings,
    scores: RsF64,
    essential: RsStrings,
    nonessential: RsStrings,
) -> *mut RsResult {
    let _profile_session = rsbagel_core::profiling::Session::new();
    let _profile_span = rsbagel_core::profiling::span("ffi_workflow");
    execute(|| {
        if genes.len != scores.len {
            return Err(FfiError::invalid("PR gene/score lengths differ"));
        }
        let genes = strings(genes, "genes")?;
        let scores = checked_slice(scores.data, scores.len, "scores")?;
        let essential = strings(essential, "essential genes")?;
        let nonessential = strings(nonessential, "nonessential genes")?;
        let result = rsbagel_core::run_pr(&genes, scores, &essential, &nonessential)?;
        let labels: Vec<&str> = result.iter().map(|row| row.gene.as_str()).collect();
        let values = result
            .iter()
            .flat_map(|row| [row.score, row.recall, row.precision, row.fdr])
            .collect();
        success(
            json!({"workflow":"pr","genes":labels,"columns":["score","recall","precision","fdr"]}),
            vec![buffer("precision_recall", result.len(), 4, values)?],
        )
    })
}

/// # Safety
/// A nonnull handle must have been returned by this library and remain live.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_result_status(result: *const RsResult) -> i32 {
    if result.is_null() {
        RSB_INVALID_ARGUMENT
    } else {
        (*result).status
    }
}
/// # Safety
/// Same live-handle requirement as `rsbagel_result_status`; returned bytes are
/// borrowed until result_free and are UTF-8, length-delimited, not NUL-terminated.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_result_error(result: *const RsResult) -> RsBytes {
    if result.is_null() {
        return RsBytes {
            data: ptr::null(),
            len: 0,
        };
    }
    let error = &(*result).error;
    RsBytes {
        data: if error.is_empty() {
            ptr::null()
        } else {
            error.as_ptr()
        },
        len: error.len(),
    }
}
/// # Safety
/// Same ownership/lifetime requirements as `rsbagel_result_error`.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_result_json(result: *const RsResult) -> RsBytes {
    if result.is_null() {
        return RsBytes {
            data: ptr::null(),
            len: 0,
        };
    }
    let metadata = &(*result).metadata;
    RsBytes {
        data: if metadata.is_empty() {
            ptr::null()
        } else {
            metadata.as_ptr()
        },
        len: metadata.len(),
    }
}
/// # Safety
/// A nonnull handle must be a live library result.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_result_buffer_count(result: *const RsResult) -> usize {
    if result.is_null() {
        0
    } else {
        (*result).buffers.len()
    }
}
/// # Safety
/// The handle remains live until all borrowed numeric views are no longer used.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_result_buffer(
    result: *const RsResult,
    index: usize,
) -> RsNumericView {
    let empty = RsNumericView {
        data: ptr::null(),
        len: 0,
        rows: 0,
        columns: 0,
    };
    if result.is_null() {
        return empty;
    }
    let buffers = &(*result).buffers;
    match buffers.get(index) {
        None => empty,
        Some(buffer) => RsNumericView {
            data: if buffer.values.is_empty() {
                ptr::null()
            } else {
                buffer.values.as_ptr()
            },
            len: buffer.values.len(),
            rows: buffer.rows,
            columns: buffer.columns,
        },
    }
}
/// # Safety
/// A nonnull result must be an owned library handle, freed exactly once. No
/// accessors/borrowed views may be in use concurrently with this call.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_result_free(result: *mut RsResult) {
    if !result.is_null() {
        drop(Box::from_raw(result));
    }
}
#[no_mangle]
pub extern "C" fn rsbagel_cancel_new() -> *mut RsCancelToken {
    Box::into_raw(Box::new(RsCancelToken {
        cancelled: AtomicBool::new(false),
    }))
}
/// # Safety
/// A nonnull token must remain a live library-owned token. Request is thread-safe.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_cancel_request(token: *const RsCancelToken) {
    if !token.is_null() {
        (*token).cancelled.store(true, Ordering::Relaxed);
    }
}
/// # Safety
/// Reset a live token only before starting the next run, not during an active run.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_cancel_reset(token: *const RsCancelToken) {
    if !token.is_null() {
        (*token).cancelled.store(false, Ordering::Relaxed);
    }
}
/// # Safety
/// Free a library token exactly once after all active calls/requests complete.
#[no_mangle]
pub unsafe extern "C" fn rsbagel_cancel_free(token: *mut RsCancelToken) {
    if !token.is_null() {
        drop(Box::from_raw(token));
    }
}
