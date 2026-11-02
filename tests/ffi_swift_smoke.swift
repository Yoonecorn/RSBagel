import Darwin

final class OwnedText {
    let pointer: UnsafeMutablePointer<CChar>
    let length: Int
    init(_ text: String) { pointer = strdup(text)!; length = text.utf8.count }
    deinit { free(pointer) }
    var view: RsbBytes { RsbBytes(data: UnsafeRawPointer(pointer).assumingMemoryBound(to: UInt8.self), len: length) }
}
let ids = [OwnedText("g0"), OwnedText("g1")]
let genes = [OwnedText("A"), OwnedText("B")]
let columns = [OwnedText("T0"), OwnedText("R1")]
let idViews = ids.map { $0.view }
let geneViews = genes.map { $0.view }
let columnViews = columns.map { $0.view }
let values: [Double] = [100, 50, 200, 90]
let config = OwnedText("{\"controls\":[0],\"control_labels\":[\"T0\"],\"min_reads\":0,\"pseudocount\":5}")
precondition(rsbagel_abi_version() == 1)
idViews.withUnsafeBufferPointer { i in
    geneViews.withUnsafeBufferPointer { g in
        columnViews.withUnsafeBufferPointer { c in
            values.withUnsafeBufferPointer { v in
                let matrix = RsbMatrix(row_ids: RsbStrings(data: i.baseAddress, len: i.count), row_genes: RsbStrings(data: g.baseAddress, len: g.count), columns: RsbStrings(data: c.baseAddress, len: c.count), values: RsbF64(data: v.baseAddress, len: v.count))
                let result = rsbagel_fc(matrix, config.view)!
                defer { rsbagel_result_free(result) }
                precondition(rsbagel_result_status(result) == RSB_OK)
                precondition(rsbagel_result_buffer_count(result) >= 2)
                let fc = rsbagel_result_buffer(result, 0)
                precondition(fc.rows == 2 && fc.columns == 1 && fc.len == 2)
                precondition(fc.data![0].isFinite && fc.data![1].isFinite)
            }
        }
    }
}
print("Swift C-ABI in-memory FC and ownership PASS")
