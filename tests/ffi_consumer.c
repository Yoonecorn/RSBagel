/* Actual C client: no files or Python subprocesses; link to librsbagel_ffi. */
#include "rsbagel.h"
#include <assert.h>
#include <math.h>
#include <stdatomic.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

enum { GENES = 120, GUIDES_PER_GENE = 3, ROWS = GENES * GUIDES_PER_GENE };
static RsbBytes bytes(const char *text) {
    RsbBytes result = {(const uint8_t *)text, strlen(text)};
    return result;
}
static RsbStrings strings(const RsbBytes *data, size_t length) {
    RsbStrings result = {data, length}; return result;
}
static void require_status(RsbResult *result, int expected) {
    assert(result != NULL);
    int status = rsbagel_result_status(result);
    if (status != expected) {
        RsbBytes error = rsbagel_result_error(result);
        fprintf(stderr, "unexpected status %d (expected %d): %.*s\n", status,
                expected, (int)error.len, error.data ? (const char *)error.data : "");
        abort();
    }
}
typedef struct { atomic_size_t calls; atomic_size_t total; } ProgressState;
static int32_t observe(void *context, size_t completed, size_t total) {
    ProgressState *state = context;
    assert(completed <= total || total == 0);
    atomic_fetch_add(&state->calls, 1);
    atomic_store(&state->total, total);
    return 0;
}
static int32_t stop(void *context, size_t completed, size_t total) {
    (void)context; (void)completed; (void)total; return 1;
}
int main(void) {
    assert(rsbagel_abi_version() == 1);
    assert(rsbagel_result_status(NULL) == RSB_INVALID_ARGUMENT);
    assert(rsbagel_result_buffer(NULL, 0).len == 0);
    assert(rsbagel_result_json(NULL).len == 0);
    rsbagel_result_free(NULL); rsbagel_cancel_request(NULL);
    rsbagel_cancel_reset(NULL); rsbagel_cancel_free(NULL);

    char gene_text[GENES][16], guide_text[ROWS][24];
    RsbBytes genes[GENES], row_ids[ROWS], row_genes[ROWS];
    RsbBytes counts_columns[] = {{0, 0}, {0, 0}, {0, 0}};
    counts_columns[0] = bytes("T0"); counts_columns[1] = bytes("R1");
    counts_columns[2] = bytes("R2");
    double counts[ROWS * 3];
    for (size_t gene = 0; gene < GENES; ++gene) {
        snprintf(gene_text[gene], sizeof gene_text[gene], "G%03zu", gene);
        genes[gene] = bytes(gene_text[gene]);
        for (size_t guide = 0; guide < GUIDES_PER_GENE; ++guide) {
            size_t row = gene * GUIDES_PER_GENE + guide;
            snprintf(guide_text[row], sizeof guide_text[row], "G%03zu_g%zu", gene, guide);
            row_ids[row] = bytes(guide_text[row]); row_genes[row] = genes[gene];
            double base = 180 + (gene * 37 + guide * 71) % 650;
            double depletion = gene < 40 ? -2.3 : gene < 80 ? 0.15 : -0.6;
            double deviation = ((int)((gene * 19 + guide * 13) % 47) - 23) / 30.0;
            counts[row * 3] = base;
            counts[row * 3 + 1] = fmax(1, nearbyint(base * pow(2, depletion + deviation)));
            counts[row * 3 + 2] = fmax(1, nearbyint(base * pow(2, depletion + deviation + ((int)((gene + guide * 3) % 9) - 4) / 13.0)));
        }
    }
    RsbMatrix input = {strings(row_ids, ROWS), strings(row_genes, ROWS),
                       strings(counts_columns, 3), {counts, ROWS * 3}};
    RsbBytes no_config = {NULL, 0};
    RsbResult *fc = rsbagel_fc(input, no_config);
    require_status(fc, RSB_OK);
    assert(rsbagel_result_buffer_count(fc) == 3);
    RsbNumericView fold_change = rsbagel_result_buffer(fc, 0);
    RsbNumericView normalized = rsbagel_result_buffer(fc, 1);
    RsbNumericView cli_bf_input = rsbagel_result_buffer(fc, 2);
    assert(cli_bf_input.rows == ROWS && cli_bf_input.columns == 2);
    assert(fold_change.rows == ROWS && fold_change.columns == 2);
    assert(normalized.rows == ROWS && normalized.columns == 3);
    assert(isfinite(fold_change.data[0]) && isfinite(normalized.data[0]));
    assert(rsbagel_result_buffer(fc, 99).len == 0);
    assert(rsbagel_result_json(fc).len > 0);

    RsbMatrix bad = input;
    bad.values.len = 1;
    RsbResult *error = rsbagel_fc(bad, no_config);
    require_status(error, RSB_INVALID_ARGUMENT); rsbagel_result_free(error);
    bad = input; bad.row_ids.data = NULL;
    error = rsbagel_fc(bad, no_config);
    require_status(error, RSB_INVALID_ARGUMENT); rsbagel_result_free(error);
    bad = input; bad.values.data = (const double *)((const uint8_t *)counts + 1);
    error = rsbagel_fc(bad, no_config);
    require_status(error, RSB_INVALID_ARGUMENT); rsbagel_result_free(error);
    bad = input; bad.row_ids.len = SIZE_MAX; bad.row_genes.len = SIZE_MAX;
    error = rsbagel_fc(bad, no_config);
    require_status(error, RSB_INVALID_ARGUMENT); rsbagel_result_free(error);
    const uint8_t invalid_utf8[] = {0xff};
    const RsbBytes invalid_id[] = {{invalid_utf8, 1}};
    bad = input; bad.row_ids = strings(invalid_id, 1);
    bad.row_genes.len = 1; bad.values.len = 3;
    error = rsbagel_fc(bad, no_config);
    require_status(error, RSB_INVALID_ARGUMENT); rsbagel_result_free(error);
    error = rsbagel_fc(input, bytes("{bad"));
    require_status(error, RSB_INVALID_ARGUMENT); rsbagel_result_free(error);

    /* NaN and infinity stay in native f64 buffers instead of JSON conversion. */
    double zero_counts[] = {0, 0, 0, 0};
    RsbMatrix nonfinite = {strings(row_ids, 2), strings(row_genes, 2),
                           strings(counts_columns, 2), {zero_counts, 4}};
    RsbResult *nonfinite_result = rsbagel_fc(nonfinite, bytes("{\"pseudocount\":0}"));
    require_status(nonfinite_result, RSB_OK);
    assert(isnan(rsbagel_result_buffer(nonfinite_result, 0).data[0]));
    rsbagel_result_free(nonfinite_result);
    double infinite_counts[] = {0, 10, 10, 10};
    nonfinite.values.data = infinite_counts;
    nonfinite_result = rsbagel_fc(nonfinite, bytes("{\"pseudocount\":0}"));
    require_status(nonfinite_result, RSB_OK);
    assert(isinf(rsbagel_result_buffer(nonfinite_result, 0).data[0]));
    assert(rsbagel_result_buffer(nonfinite_result, 0).data[0] > 0);
    rsbagel_result_free(nonfinite_result);

    RsbBytes fc_columns[] = {{0,0},{0,0}};
    fc_columns[0] = bytes("R1"); fc_columns[1] = bytes("R2");
    RsbMatrix fc_input = {strings(row_ids, ROWS), strings(row_genes, ROWS),
                         strings(fc_columns, 2), {cli_bf_input.data, cli_bf_input.len}};
    RsbBytes config = bytes("{\"columns\":[0,1],\"column_labels\":[\"R1\",\"R2\"],\"no_resampling\":true,\"rng_mode\":\"numpy\",\"seed\":42,\"threads\":1}");
    RsbStrings essential = strings(genes, 40), nonessential = strings(genes + 40, 40);
    /* Omitted BF config uses production NumPy sampling. Keep resampling enabled
       so this checks the RNG-dependent workflow, not only a default field. */
    RsbResult *default_bf = rsbagel_bf(fc_input, essential, nonessential, no_config,
                                     NULL, NULL, NULL);
    RsbResult *numpy_bf = rsbagel_bf(fc_input, essential, nonessential,
                                   bytes("{\"rng_mode\":\"numpy\"}"), NULL, NULL, NULL);
    require_status(default_bf, RSB_OK); require_status(numpy_bf, RSB_OK);
    assert(rsbagel_result_buffer_count(default_bf) == rsbagel_result_buffer_count(numpy_bf));
    for (size_t index = 0; index < rsbagel_result_buffer_count(default_bf); ++index) {
        RsbNumericView actual = rsbagel_result_buffer(default_bf, index);
        RsbNumericView expected = rsbagel_result_buffer(numpy_bf, index);
        assert(actual.rows == expected.rows && actual.columns == expected.columns);
        assert(actual.len == expected.len);
        assert(memcmp(actual.data, expected.data, actual.len * sizeof(double)) == 0);
    }
    RsbBytes default_json = rsbagel_result_json(default_bf);
    RsbBytes numpy_json = rsbagel_result_json(numpy_bf);
    assert(default_json.len == numpy_json.len);
    assert(memcmp(default_json.data, numpy_json.data, default_json.len) == 0);
    rsbagel_result_free(default_bf); rsbagel_result_free(numpy_bf);
    ProgressState progress = {ATOMIC_VAR_INIT(0), ATOMIC_VAR_INIT(0)};
    RsbResult *bf = rsbagel_bf(fc_input, essential, nonessential, config,
                              observe, &progress, NULL);
    require_status(bf, RSB_OK);
    assert(atomic_load(&progress.calls) >= 2);
    assert(atomic_load(&progress.total) == 1);
    assert(rsbagel_result_buffer_count(bf) == 5);
    RsbNumericView gene_scores = rsbagel_result_buffer(bf, 0);
    assert(gene_scores.rows == GENES && gene_scores.columns == 7);
    RsbNumericView cli_pr_input = rsbagel_result_buffer(bf, 4);
    assert(cli_pr_input.rows == GENES && cli_pr_input.columns == 1);
    double scores[GENES];
    for (size_t gene = 0; gene < GENES; ++gene) {
        scores[gene] = cli_pr_input.data[gene];
        assert(isfinite(scores[gene]));
    }
    RsbResult *pr = rsbagel_pr(strings(genes, GENES), (RsbF64){scores, GENES},
                              essential, nonessential);
    require_status(pr, RSB_OK);
    RsbNumericView precision_recall = rsbagel_result_buffer(pr, 0);
    assert(precision_recall.rows == GENES && precision_recall.columns == 4);
    for (size_t row = 1; row < GENES; ++row)
        assert(precision_recall.data[(row - 1) * 4] >= precision_recall.data[row * 4]);
    assert(precision_recall.data[(GENES - 1) * 4 + 1] == 1.0);

    RsbCancelToken *token = rsbagel_cancel_new();
    assert(token != NULL); rsbagel_cancel_request(token);
    RsbResult *cancelled = rsbagel_bf(fc_input, essential, nonessential, config,
                                    NULL, NULL, token);
    require_status(cancelled, RSB_CANCELLED); rsbagel_result_free(cancelled);
    rsbagel_cancel_reset(token);
    cancelled = rsbagel_bf(fc_input, essential, nonessential, config, stop, NULL, token);
    require_status(cancelled, RSB_CANCELLED); rsbagel_result_free(cancelled);
    rsbagel_cancel_free(token);

    rsbagel_result_free(pr); rsbagel_result_free(bf); rsbagel_result_free(fc);
    puts("ffi_consumer: fc/bf/pr, ownership, pointer validation, native nonfinites, progress and cancellation PASS");
    return 0;
}
