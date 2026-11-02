#ifndef RSBAGEL_H
#define RSBAGEL_H
#include <stddef.h>
#include <stdint.h>
#ifdef __cplusplus
extern "C" {
#endif
/* ABI version 1. Strings are UTF-8 pointer/length pairs, not NUL terminated.
 * Zero length permits a NULL data pointer. Nonempty inputs must be aligned,
 * valid allocations of the declared length, live until the call returns.
 * Inputs are borrowed synchronously and copied into the in-memory Rust model.
 */
typedef struct { const uint8_t *data; size_t len; } RsbBytes;
typedef struct { const RsbBytes *data; size_t len; } RsbStrings;
typedef struct { const double *data; size_t len; } RsbF64;
typedef struct {
    RsbStrings row_ids;
    RsbStrings row_genes;
    RsbStrings columns;
    RsbF64 values; /* row-major; row_ids.len * columns.len elements */
} RsbMatrix;
typedef struct { const double *data; size_t len, rows, columns; } RsbNumericView;
typedef struct RsResult RsbResult;
typedef struct RsCancelToken RsbCancelToken;
enum {
    RSB_OK = 0, RSB_INVALID_ARGUMENT = 1, RSB_NUMERICAL_ERROR = 2,
    RSB_LEGACY_FAILURE = 3, RSB_CANCELLED = 4, RSB_IO_ERROR = 5,
    RSB_INTERNAL_ERROR = 6
};
/* Callback may run on a worker thread. Context/callback must be thread-safe and
 * must not unwind or free active inputs/tokens. Return nonzero to cancel.
 * The initial callback uses (0, 0) before the core reports its iteration total.
 */
typedef int32_t (*RsbProgress)(void *context, size_t completed, size_t total);
uint32_t rsbagel_abi_version(void);
/* Empty config bytes select defaults. FC keys: controls (zero-based numeric
 * columns), control_labels, min_reads, pseudocount. BF follows BfConfig JSON.
 * Every returned handle, including an error handle, must be freed exactly once.
 * There is no disk I/O or subprocess execution in these functions.
 */
/* FC buffers0/1 remain raw FC/normalized reads. Optional buffer2 is
 * cli_bf_input, rounded exactly at the upstream 3-decimal FC writer boundary.
 * Check cli_bf_input_available in JSON; invalid empty-field transitions are
 * flagged without inventing a successful BF input. */
RsbResult *rsbagel_fc(RsbMatrix input, RsbBytes config_json);
/* BF raw buffers0..3 are unchanged. Optional buffer4 cli_pr_input applies
 * the source BF-writer scaling/3-decimal/alphabetical boundary. Its gene
 * identities are cli_pr_genes in JSON, not necessarily raw-buffer row order.
 * Use ready inputs directly; do not multiply their replicate factor again. */
RsbResult *rsbagel_bf(RsbMatrix input, RsbStrings essential, RsbStrings nonessential,
    RsbBytes config_json, RsbProgress callback, void *context,
    const RsbCancelToken *cancel_token);
RsbResult *rsbagel_pr(RsbStrings genes, RsbF64 scores, RsbStrings essential,
    RsbStrings nonessential);
/* A nonnull result must be a live handle from this library. Accessors borrow;
 * they neither allocate nor transfer ownership. Views expire at result_free.
 * NULL result: INVALID_ARGUMENT status, empty views and zero buffer count.
 * Out-of-range buffer index: empty view. Do not free borrowed view pointers.
 * Numeric views preserve IEEE-754 NaN/Inf. JSON may encode nonfinite as null;
 * metadata records optional-field presence separately from numeric contents.
 */
int32_t rsbagel_result_status(const RsbResult *result);
RsbBytes rsbagel_result_error(const RsbResult *result);
RsbBytes rsbagel_result_json(const RsbResult *result);
size_t rsbagel_result_buffer_count(const RsbResult *result);
RsbNumericView rsbagel_result_buffer(const RsbResult *result, size_t index);
void rsbagel_result_free(RsbResult *result); /* NULL is a no-op */
/* Cancellation token must outlive each run using it. Request is concurrent-safe;
 * reset only before a new run, and free after all active runs/requests finish.
 */
RsbCancelToken *rsbagel_cancel_new(void);
void rsbagel_cancel_request(const RsbCancelToken *token);
void rsbagel_cancel_reset(const RsbCancelToken *token);
void rsbagel_cancel_free(RsbCancelToken *token); /* NULL is a no-op */
#ifdef __cplusplus
}
#endif
#endif
