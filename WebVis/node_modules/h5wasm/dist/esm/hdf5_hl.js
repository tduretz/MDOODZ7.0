import ModuleFactory from './hdf5_util.js';
export var Module; //: H5WasmModule = null;
export var FS = null;
const ready = ModuleFactory({ noInitialRun: true }).then(result => { Module = result; FS = Module.FS; return Module; });
export { ready };
export const ACCESS_MODES = {
    "r": "H5F_ACC_RDONLY",
    "a": "H5F_ACC_RDWR",
    "w": "H5F_ACC_TRUNC",
    "x": "H5F_ACC_EXCL",
    "Sa": "ACC_SWMR_APPEND",
    "Sr": "H5F_ACC_SWMR_READ"
};
export const LIBVER_BOUNDS_MAP = {
    "earliest": "H5F_LIBVER_EARLIEST",
    "v108": "H5F_LIBVER_V18",
    "v110": "H5F_LIBVER_V110",
    "v112": "H5F_LIBVER_V112",
    "v114": "H5F_LIBVER_V114",
    "v200": "H5F_LIBVER_V200",
    "latest": "H5F_LIBVER_LATEST"
};
function normalizePath(path) {
    if (path == "/") {
        return path;
    }
    // replace multiple path separators with single
    path = path.replace(/\/(\/)+/g, '/');
    // strip end slashes
    path = path.replace(/(\/)+$/, '');
    return path;
}
function dirname(path) {
    // adapted from dirname function in posixpath.py
    const sep = "/";
    const sep_index = path.lastIndexOf(sep) + 1;
    let head = path.slice(0, sep_index);
    if (head && head !== sep.repeat(head.length)) {
        // strip end slashes
        head = head.replace(/(\/)+$/, '');
    }
    return head;
}
export function convertToLibverString(value) {
    // Create reverse mapping from constant name to key
    for (const [key, constantName] of Object.entries(LIBVER_BOUNDS_MAP)) {
        const moduleValue = Module[constantName];
        if (moduleValue === value) {
            return key;
        }
    }
    throw new Error(`No matching libver string found for value: ${value}`);
}
function parseLibverString(s) {
    const normalized = s.trim().toLowerCase();
    const constant_name = LIBVER_BOUNDS_MAP[normalized];
    if (!constant_name) {
        throw new Error(`Invalid libver value: "${s}". Use: ${Object.keys(LIBVER_BOUNDS_MAP).join(', ')}`);
    }
    const value = Module[constant_name];
    if (value === undefined) {
        throw new Error(`Module constant ${constant_name} is not defined`);
    }
    return value;
}
function check_malloc(nbytes) {
    const max_memory = Module.MAXIMUM_MEMORY;
    if (nbytes > max_memory) {
        throw new Error(`Requested allocation of ${nbytes} bytes exceeds maximum memory of ${max_memory} bytes`);
    }
    const safe_nbytes = Number(nbytes);
    const ptr = Module._malloc(safe_nbytes);
    if (ptr === 0) {
        throw new Error(`Memory allocation of ${safe_nbytes} bytes failed`);
    }
    return ptr;
}
function get_attr(file_id, obj_name, attr_name, json_compatible = false) {
    let metadata = Module.get_attribute_metadata(file_id, obj_name, attr_name);
    if (!metadata.shape) {
        return null;
    }
    const nbytes = metadata.size * metadata.total_size;
    let data_ptr = check_malloc(nbytes);
    var processed;
    try {
        Module.get_attribute_data(file_id, obj_name, attr_name, BigInt(data_ptr));
        let data = Module.HEAPU8.slice(data_ptr, data_ptr + nbytes);
        processed = process_data(data, metadata, json_compatible);
    }
    finally {
        if (metadata.vlen) {
            Module.reclaim_vlen_memory(file_id, obj_name, attr_name, BigInt(data_ptr));
        }
        Module._free(data_ptr);
    }
    if (json_compatible) {
        return processed;
    }
    return processed;
}
function getAccessor(type, size, signed) {
    if (type === 0) {
        if (size === 8) {
            return (signed) ? BigInt64Array : BigUint64Array;
        }
        else if (size === 4) {
            return (signed) ? Int32Array : Uint32Array;
        }
        else if (size === 2) {
            return (signed) ? Int16Array : Uint16Array;
        }
        else { // size === 1
            return (signed) ? Int8Array : Uint8Array;
        }
    }
    else { // type ==== 1 (float)
        if (size === 8) {
            return Float64Array;
        }
        else if (size === 4) {
            return Float32Array;
        }
        else {
            throw new Error(`Float${size * 8} not supported`);
        }
    }
}
function process_data(data, metadata, json_compatible = false) {
    // (for data coming out of Module)
    // If an appropriate TypedArray container can be constructed, it will
    // but otherwise returns Uint8Array raw bytes as loaded.
    let output_data;
    let { shape, type } = metadata;
    let known_type = true;
    // let length: number;
    if (type === Module.H5T_class_t.H5T_STRING.value) {
        if (metadata.vlen) {
            let output = [];
            let reader = (metadata.cset == 1) ? Module.UTF8ToString : Module.AsciiToString;
            let ptrs = new Uint32Array(data.buffer);
            for (let ptr of ptrs) {
                output.push(reader(ptr));
            }
            output_data = output;
            // length = output_data.length;
        }
        else {
            let encoding = (metadata.cset == 1) ? 'utf-8' : 'ascii';
            let decoder = new TextDecoder(encoding);
            let size = metadata.size;
            let n = Math.floor(data.byteLength / size);
            let output = [];
            for (let i = 0; i < n; i++) {
                let bytes = data.slice(i * size, (i + 1) * size);
                // truncate at first null
                const zero_match = bytes.findIndex((b) => (b === 0));
                if (zero_match > -1) {
                    bytes = bytes.slice(0, zero_match);
                }
                output.push(decoder.decode(bytes));
            }
            output_data = output;
            // length = output_data.length;
        }
    }
    else if (type === Module.H5T_class_t.H5T_INTEGER.value || type === Module.H5T_class_t.H5T_FLOAT.value) {
        const { size, signed } = metadata;
        const accessor = getAccessor(type, size, signed);
        output_data = new accessor(data.buffer);
        if (json_compatible) {
            output_data = [...output_data];
            if (accessor === BigInt64Array || accessor === BigUint64Array) {
                output_data = output_data.map(Number);
            }
        }
    }
    else if (type === Module.H5T_class_t.H5T_COMPOUND.value) {
        const { size, compound_type } = metadata;
        let n = Math.floor(data.byteLength / size);
        let output = [];
        for (let i = 0; i < n; i++) {
            let row = [];
            let row_data = data.slice(i * size, (i + 1) * size);
            for (let member of compound_type.members) {
                let member_data = row_data.slice(member.offset, member.offset + member.size);
                row.push(process_data(member_data, member, json_compatible));
            }
            output.push(row);
        }
        output_data = output;
    }
    else if (type === Module.H5T_class_t.H5T_ARRAY.value) {
        const { array_type } = metadata;
        shape = shape.concat(array_type.shape);
        array_type.shape = shape;
        // always convert ARRAY types to base JS types:
        output_data = process_data(data, array_type, true);
        if (isIterable(output_data) && typeof output_data !== "string") {
            output_data = create_nested_array(output_data, array_type.shape);
        }
    }
    else if (type === Module.H5T_class_t.H5T_ENUM.value) {
        const base_metadata = { ...metadata };
        base_metadata.type = base_metadata.enum_type.type;
        output_data = process_data(data, base_metadata, json_compatible);
        // Following the convention of h5py, treat all enum datasets where the
        // enum members are ["FALSE", "TRUE"] as boolean arrays
        if (json_compatible && isH5PYBooleanEnum(metadata.enum_type)) {
            if (isIterable(output_data)) {
                output_data = [...output_data].map((x) => !!x);
            }
            else {
                output_data = !!output_data;
            }
        }
    }
    else if (type === Module.H5T_class_t.H5T_REFERENCE.value) {
        const { ref_type, size } = metadata; // as { ref_type: 'object' | 'region', size: number };
        const cls = (ref_type === 'object') ? Reference : RegionReference;
        output_data = Array.from({ length: metadata.total_size }).map((_, i) => {
            const ref_data = data.slice(i * size, (i + 1) * size);
            return new cls(ref_data);
        });
        return output_data;
    }
    else if (type === Module.H5T_class_t.H5T_VLEN.value) {
        // Data buffer holds HDF5 `hvl_t` struct
        // https://docs.hdfgroup.org/hdf5/v1_12/structhvl__t.html
        const ptr_pairs = new Uint32Array(data.buffer); // both `size_t` length and pointer are 4-byte-long
        const ptr_pairs_length = ptr_pairs.length;
        const vlen_type = metadata.vlen_type;
        const { size } = vlen_type;
        let output = [];
        for (let p = 0; p < ptr_pairs_length; p += 2) {
            const length = ptr_pairs[p]; // number of elements in this vlen array
            const data_ptr = ptr_pairs[p + 1]; // pointer to this vlen array's data
            // Read vlen array data from memory
            const data_nbytes = length * size;
            const data = Module.HEAPU8.slice(data_ptr, data_ptr + data_nbytes);
            // Process this vlen array's data according to base datatype
            output.push(process_data(data, { ...vlen_type, shape: [length], total_size: length }, json_compatible));
        }
        output_data = output;
    }
    else {
        known_type = false;
        output_data = data;
    }
    // if metadata.shape.length == 0 or metadata.shape is undefined...
    if (known_type && (Array.isArray(output_data) || ArrayBuffer.isView(output_data)) && !shape?.length) {
        output_data = output_data[0];
    }
    if (json_compatible) {
        return output_data;
    }
    return output_data;
}
function isIterable(x) {
    return typeof x === 'object' && x !== null && Symbol.iterator in x;
}
function isH5PYBooleanEnum(enum_type) {
    return Object.keys(enum_type.members).length === 2 &&
        enum_type.members["FALSE"] === 0 &&
        enum_type.members["TRUE"] === 1;
}
const H5S_UNLIMITED = 18446744073709551615n; // not exportable by emscripten_bindings because it's outside int32 range
function prepare_data(data, metadata, shape, maxshape) {
    // for data being sent to Module
    let final_shape;
    let final_maxshape;
    if (shape === undefined || shape === null) {
        if (data instanceof Map) {
            // For maps, try to grab the length of the first column to guess the shape
            let guessed_len = 1;
            for (const col of data.values()) {
                if (col && col.length !== undefined && !(typeof col === 'string')) {
                    guessed_len = col.length;
                    break;
                }
            }
            final_shape = [BigInt(guessed_len)];
        }
        else if (data != null && data.length != null && !(typeof data === 'string')) {
            final_shape = [BigInt(data.length)];
        }
        else {
            final_shape = [];
        }
    }
    else {
        final_shape = shape.map(BigInt);
    }
    if (maxshape === undefined || maxshape === null) {
        final_maxshape = final_shape;
    }
    else {
        final_maxshape = maxshape.map((dim) => ((dim === null) ? H5S_UNLIMITED : BigInt(dim)));
    }
    data = (Array.isArray(data) || ArrayBuffer.isView(data) || data instanceof Map) ? data : [data];
    let total_size = Number(final_shape.reduce((previous, current) => current * previous, 1n));
    if (!(data instanceof Map)) {
        if (data.length != total_size) {
            throw new Error(`Error: shape ${final_shape} does not match number of elements in data`);
        }
    }
    let output;
    if (metadata.type === Module.H5T_class_t.H5T_STRING.value) {
        if (!metadata.vlen) {
            output = new Uint8Array(total_size * metadata.size);
            let encoder = new TextEncoder();
            output.fill(0);
            let offset = 0;
            for (let s of data) {
                let encoded = encoder.encode(s);
                output.set(encoded.slice(0, metadata.size), offset);
                offset += metadata.size;
            }
        }
        else {
            output = data;
        }
    }
    else if (metadata.type === Module.H5T_class_t.H5T_INTEGER.value || metadata.type === Module.H5T_class_t.H5T_FLOAT.value) {
        const { type, size, signed } = metadata;
        const accessor = getAccessor(type, size, signed);
        let typed_array;
        if (data instanceof accessor) {
            typed_array = data;
        }
        else {
            // convert...
            if (metadata.size > 4 && metadata.type === Module.H5T_class_t.H5T_INTEGER.value) {
                data = data.map(BigInt);
            }
            typed_array = new accessor(data);
        }
        output = new Uint8Array(typed_array.buffer, typed_array.byteOffset, typed_array.byteLength);
    }
    else if (metadata.type === Module.H5T_class_t.H5T_REFERENCE.value) {
        output = new Uint8Array(metadata.size * total_size);
        data.forEach((r, i) => output.set(r.ref_data, i * metadata.size));
    }
    else if (metadata.type === Module.H5T_class_t.H5T_COMPOUND.value) {
        const { size, compound_type } = metadata;
        output = new Uint8Array(total_size * size);
        const member_buffers = new Map();
        const map_data = data;
        // Recursively convert each column into its raw flat bytes
        for (const member of compound_type.members) {
            if (member.vlen) {
                throw new Error(`Writing VLEN strings inside compound types requires pointers in Wasm heap and is not currently supported.`);
            }
            const column_data = map_data.get(member.name);
            // The recursive call will validate the column length matches the total shape
            const prepared = prepare_data(column_data, member, shape, maxshape);
            member_buffers.set(member.name, prepared.data);
        }
        // Interleave the columns into the final Array of Structures buffer
        for (let i = 0; i < total_size; i++) {
            const row_offset = i * size;
            for (const member of compound_type.members) {
                const col_buf = member_buffers.get(member.name);
                const m_size = member.size;
                // Grab the bytes for the i-th element of this specific member
                const element_bytes = col_buf.subarray(i * m_size, (i + 1) * m_size);
                // Write it into the compound struct memory space
                output.set(element_bytes, row_offset + member.offset);
            }
        }
    }
    else {
        throw new Error(`data with type ${metadata.type} can not be prepared for write`);
    }
    return { data: output, shape: final_shape, maxshape: final_maxshape };
}
function map_reverse(map) {
    return new Map(Array.from(map.entries()).map(([k, v]) => [v, k]));
}
const int_fmts = new Map([[1, 'b'], [2, 'h'], [4, 'i'], [8, 'q']]);
const float_fmts = new Map([[2, 'e'], [4, 'f'], [8, 'd']]);
const fmts_float = map_reverse(float_fmts);
const fmts_int = map_reverse(int_fmts);
function metadata_to_dtype(metadata) {
    const { type, size, littleEndian, signed, compound_type, array_type, vlen } = metadata;
    if (type == Module.H5T_class_t.H5T_STRING.value) {
        let length_str = vlen ? "" : String(size);
        return `S${length_str}`;
    }
    else if (type == Module.H5T_class_t.H5T_INTEGER.value) {
        let fmt = int_fmts.get(size);
        if (fmt === undefined) {
            throw new Error(`int of size ${size} unsupported`);
        }
        if (!signed) {
            fmt = fmt.toUpperCase();
        }
        return ((littleEndian) ? "<" : ">") + fmt;
    }
    else if (type == Module.H5T_class_t.H5T_FLOAT.value) {
        let fmt = float_fmts.get(size);
        return ((littleEndian) ? "<" : ">") + fmt;
    }
    else if (type == Module.H5T_class_t.H5T_COMPOUND.value) {
        const ct = compound_type;
        return ct.members.map((member) => {
            if (member.type === Module.H5T_class_t.H5T_ARRAY.value && member.array_type) {
                return [member.name, metadata_to_dtype(member.array_type), member.array_type.shape];
            }
            return [member.name, metadata_to_dtype(member)];
        });
    }
    else if (type === Module.H5T_class_t.H5T_ARRAY.value) {
        const at = array_type;
        return [metadata_to_dtype(at), at.shape];
    }
    else if (type === Module.H5T_class_t.H5T_REFERENCE.value) {
        return (metadata.ref_type === 'object') ? "Reference" : "RegionReference";
    }
    else {
        return "unknown";
    }
}
export function dtype_to_metadata(dtype) {
    let metadata = { vlen: false, signed: false, littleEndian: true };
    if (typeof dtype === 'string') {
        // Simple string dtype: '<i8', '<f4', 'S10', 'Reference', etc.
        if (dtype === "Reference" || dtype === "RegionReference") {
            metadata.type = Module.H5T_class_t.H5T_REFERENCE.value;
            metadata.size = (dtype === "Reference") ? Module.SIZEOF_OBJ_REF : Module.SIZEOF_DSET_REGION_REF;
        }
        else {
            const match = dtype.match(/^([<>|]?)([bhiqefdsBHIQS])([0-9]*)$/);
            if (match == null) {
                throw dtype + " is not a recognized dtype";
            }
            const [, endianness, typestr, length] = match;
            metadata.littleEndian = (endianness != '>');
            if (fmts_int.has(typestr.toLowerCase())) {
                metadata.type = Module.H5T_class_t.H5T_INTEGER.value;
                metadata.size = fmts_int.get(typestr.toLowerCase());
                metadata.signed = (typestr.toLowerCase() == typestr);
            }
            else if (fmts_float.has(typestr)) {
                metadata.type = Module.H5T_class_t.H5T_FLOAT.value;
                metadata.size = fmts_float.get(typestr);
            }
            else if (typestr.toUpperCase() === 'S') {
                metadata.type = Module.H5T_class_t.H5T_STRING.value;
                metadata.size = (length == "") ? 4 : parseInt(length, 10);
                metadata.vlen = (length == "");
            }
            else {
                throw "should never happen";
            }
        }
    }
    else if (Array.isArray(dtype)) {
        if (Array.isArray(dtype[0])) {
            // CompoundDtype: [[name, dtype], [name, dtype, shape], ...]
            // Packed layout (no padding): offsets are a simple running sum of sizes.
            // HDF5 compound types are explicitly laid out, so packing is valid and
            // matches what h5py produces by default.
            const compound = dtype;
            let offset = 0;
            const members = compound.map(([name, memberDtype, arrayShape]) => {
                let memberMeta;
                if (arrayShape !== undefined) {
                    const baseMeta = dtype_to_metadata(memberDtype);
                    const nelems = arrayShape.reduce((a, b) => a * b, 1);
                    memberMeta = {
                        ...baseMeta,
                        type: Module.H5T_class_t.H5T_ARRAY.value,
                        array_type: baseMeta,
                        shape: arrayShape,
                        size: baseMeta.size * nelems,
                    };
                }
                else {
                    memberMeta = dtype_to_metadata(memberDtype);
                }
                const member = { ...memberMeta, name, offset };
                offset += memberMeta.size;
                return member;
            });
            metadata.type = Module.H5T_class_t.H5T_COMPOUND.value;
            metadata.compound_type = { members, nmembers: members.length };
            metadata.size = offset;
        }
        else {
            // ArrayDtype: [dtype_str, shape]
            const [base_dtype_str, shape] = dtype;
            const baseMeta = dtype_to_metadata(base_dtype_str);
            const nelems = shape.reduce((a, b) => a * b, 1);
            metadata.type = Module.H5T_class_t.H5T_ARRAY.value;
            metadata.array_type = baseMeta;
            metadata.shape = shape;
            metadata.size = baseMeta.size * nelems;
        }
    }
    else {
        // Legacy object forms: { compound_type: ... } | { array_type: ... }
        if ('compound_type' in dtype) {
            const ct = dtype.compound_type;
            metadata.type = Module.H5T_class_t.H5T_COMPOUND.value;
            metadata.compound_type = ct;
            metadata.size = ct.members.reduce((sum, m) => Math.max(sum, m.offset + m.size), 0);
        }
        else {
            const at = dtype.array_type;
            const shape = at.shape;
            const nelems = shape.reduce((a, b) => a * b, 1);
            metadata.type = Module.H5T_class_t.H5T_ARRAY.value;
            metadata.array_type = at;
            metadata.shape = shape;
            metadata.size = at.size * nelems;
        }
    }
    return metadata;
}
const TypedArray_to_dtype = new Map([
    ['Uint8Array', '<B'],
    ['Uint8ClampedArray', '<B'],
    ['Uint16Array', '<H'],
    ['Uint32Array', '<I'],
    ['BigUint64Array', '<Q'],
    ['Int8Array', '<b'],
    ['Int16Array', '<h'],
    ['Int32Array', '<i'],
    ['BigInt64Array', '<q'],
    ['Float32Array', '<f'],
    ['Float64Array', '<d']
]);
function guess_dtype(data) {
    if (ArrayBuffer.isView(data)) {
        const dtype = TypedArray_to_dtype.get(data.constructor.name);
        if (dtype === undefined) {
            throw new Error("DataView not supported directly for write");
        }
        return dtype;
    }
    else {
        // then it is an array or a single value...
        const arr_data = ((Array.isArray(data)) ? data : [data]);
        if (arr_data.every(Number.isInteger)) {
            return '<i'; // default integer type: Int32
        }
        else if (arr_data.every((d) => (typeof d == 'number'))) {
            return '<d'; // default float type: Float64
        }
        else if (arr_data.every((d) => (typeof d == 'string'))) {
            return 'S';
        }
        else if (arr_data.every((d) => d instanceof RegionReference)) {
            return 'RegionReference';
        }
        else if (arr_data.every((d) => d instanceof Reference)) {
            return 'Reference';
        }
    }
    throw new Error("unguessable type for data");
}
export function guess_metadata(data) {
    const baseMeta = { vlen: false, signed: false, littleEndian: true };
    if (data instanceof Map) {
        let offset = 0;
        const members = [];
        for (const [name, value] of data.entries()) {
            // Recursively guess the metadata for this specific column/member
            const memberMeta = guess_metadata(value);
            members.push({
                ...memberMeta,
                name,
                offset
            });
            // Increment the running offset by the size of this member
            offset += memberMeta.size;
        }
        return {
            ...baseMeta,
            type: Module.H5T_class_t.H5T_COMPOUND.value,
            compound_type: { members, nmembers: members.length },
            size: offset // Total size of the packed compound element
        };
    }
    if (ArrayBuffer.isView(data)) {
        const cname = data.constructor.name;
        switch (cname) {
            case 'Float32Array': return { ...baseMeta, type: Module.H5T_class_t.H5T_FLOAT.value, size: 4 };
            case 'Float64Array': return { ...baseMeta, type: Module.H5T_class_t.H5T_FLOAT.value, size: 8 };
            case 'Int8Array': return { ...baseMeta, type: Module.H5T_class_t.H5T_INTEGER.value, size: 1, signed: true };
            case 'Int16Array': return { ...baseMeta, type: Module.H5T_class_t.H5T_INTEGER.value, size: 2, signed: true };
            case 'Int32Array': return { ...baseMeta, type: Module.H5T_class_t.H5T_INTEGER.value, size: 4, signed: true };
            case 'BigInt64Array': return { ...baseMeta, type: Module.H5T_class_t.H5T_INTEGER.value, size: 8, signed: true };
            case 'Uint8Array':
            case 'Uint8ClampedArray': return { ...baseMeta, type: Module.H5T_class_t.H5T_INTEGER.value, size: 1 };
            case 'Uint16Array': return { ...baseMeta, type: Module.H5T_class_t.H5T_INTEGER.value, size: 2 };
            case 'Uint32Array': return { ...baseMeta, type: Module.H5T_class_t.H5T_INTEGER.value, size: 4 };
            case 'BigUint64Array': return { ...baseMeta, type: Module.H5T_class_t.H5T_INTEGER.value, size: 8 };
            default: throw new Error("DataView not supported directly for write");
        }
    }
    else {
        // then it is an array or a single value...
        const arr_data = ((Array.isArray(data)) ? data : [data]);
        if (arr_data.every(Number.isInteger)) {
            return { ...baseMeta, type: Module.H5T_class_t.H5T_INTEGER.value, size: 4, signed: true };
        }
        else if (arr_data.every((d) => (typeof d == 'number'))) {
            return { ...baseMeta, type: Module.H5T_class_t.H5T_FLOAT.value, size: 8 };
        }
        else if (arr_data.every((d) => (typeof d == 'string'))) {
            return { ...baseMeta, type: Module.H5T_class_t.H5T_STRING.value, size: 4, vlen: true };
        }
        else if (arr_data.every((d) => d instanceof RegionReference)) {
            return { ...baseMeta, type: Module.H5T_class_t.H5T_REFERENCE.value, size: Module.SIZEOF_DSET_REGION_REF };
        }
        else if (arr_data.every((d) => d instanceof Reference)) {
            return { ...baseMeta, type: Module.H5T_class_t.H5T_REFERENCE.value, size: Module.SIZEOF_OBJ_REF };
        }
    }
    throw new Error("unguessable type for data");
}
var OBJECT_TYPE;
(function (OBJECT_TYPE) {
    OBJECT_TYPE["DATASET"] = "Dataset";
    OBJECT_TYPE["GROUP"] = "Group";
    OBJECT_TYPE["BROKEN_SOFT_LINK"] = "BrokenSoftLink";
    OBJECT_TYPE["EXTERNAL_LINK"] = "ExternalLink";
    OBJECT_TYPE["DATATYPE"] = "Datatype";
    OBJECT_TYPE["REFERENCE"] = "Reference";
    OBJECT_TYPE["REGION_REFERENCE"] = "RegionReference";
})(OBJECT_TYPE || (OBJECT_TYPE = {}));
export class BrokenSoftLink {
    constructor(target) {
        this.type = OBJECT_TYPE.BROKEN_SOFT_LINK;
        this.target = target;
    }
}
export class ExternalLink {
    constructor(filename, obj_path) {
        this.type = OBJECT_TYPE.EXTERNAL_LINK;
        this.filename = filename;
        this.obj_path = obj_path;
    }
}
export class Reference {
    constructor(ref_data) {
        this.ref_data = ref_data;
    }
}
export class RegionReference extends Reference {
}
export class Attribute {
    constructor(file_id, path, name) {
        this.file_id = file_id;
        this.path = path;
        this.name = name;
        const metadata = Module.get_attribute_metadata(file_id, path, name);
        this.metadata = metadata;
        this.dtype = metadata_to_dtype(metadata);
        this.shape = metadata.shape;
    }
    get value() {
        if (typeof this._value === "undefined") {
            this._value = get_attr(this.file_id, this.path, this.name, false);
        }
        return this._value;
    }
    get json_value() {
        if (typeof this._json_value === "undefined") {
            this._json_value = get_attr(this.file_id, this.path, this.name, true);
        }
        return this._json_value;
    }
    to_array() {
        const { json_value, metadata } = this;
        const { shape } = metadata;
        if (!isIterable(json_value) || typeof json_value === "string" || shape === null) {
            return json_value;
        }
        return create_nested_array(json_value, shape);
    }
}
class HasAttrs {
    get attrs() {
        let attr_names = Module.get_attribute_names(this.file_id, this.path);
        let attrs = {};
        const { file_id, path } = this;
        for (let name of attr_names) {
            Object.defineProperty(attrs, name, {
                get: () => (new Attribute(file_id, path, name)),
                enumerable: true
            });
        }
        return attrs;
    }
    get root() {
        return new Group(this.file_id, '/');
    }
    get parent() {
        return this.root.get(dirname(this.path));
    }
    get_attribute(name, json_compatible = false) {
        return get_attr(this.file_id, this.path, name, json_compatible);
    }
    create_attribute(name, data, shape, dtype) {
        let metadata = dtype ? dtype_to_metadata(dtype) : guess_metadata(data);
        if (!metadata.littleEndian) {
            throw new Error("create_attribute with big-endian dtype is not supported");
        }
        const { data: prepared_data, shape: guessed_shape } = prepare_data(data, metadata, shape, shape);
        const final_shape = shape?.map(BigInt) ?? guessed_shape;
        const attr_metadata = {
            ...metadata,
            shape: final_shape,
            maxshape: final_shape // setup_dataset handles maxshape; attributes can just mirror shape
        };
        if (metadata.vlen) {
            Module.create_vlen_str_attribute(this.file_id, this.path, name, prepared_data, attr_metadata);
        }
        else {
            let data_ptr = check_malloc(prepared_data.byteLength);
            try {
                Module.HEAPU8.set(prepared_data, data_ptr);
                Module.create_attribute(this.file_id, this.path, name, BigInt(data_ptr), attr_metadata);
            }
            finally {
                Module._free(data_ptr);
            }
        }
    }
    delete_attribute(name) {
        // returns non-zero value if delete failed.
        return Module.delete_attribute(this.file_id, this.path, name);
    }
    create_reference() {
        const ref_data = Module.create_object_reference(this.file_id, this.path);
        return new Reference(ref_data);
    }
    dereference(ref) {
        const is_region = (ref instanceof RegionReference);
        const name = Module.get_referenced_name(this.file_id, ref.ref_data, !is_region);
        const target = this.root.get(name);
        return (is_region) ? new DatasetRegion(target, ref) : target;
    }
}
export class Datatype extends HasAttrs {
    constructor(file_id, path) {
        super();
        this.file_id = file_id;
        this.path = path;
        this.type = OBJECT_TYPE.DATATYPE;
    }
    get metadata() {
        return Module.get_datatype_metadata(this.file_id, this.path);
    }
}
export class Group extends HasAttrs {
    constructor(file_id, path) {
        super();
        this.path = path;
        this.file_id = file_id;
        this.type = OBJECT_TYPE.GROUP;
    }
    keys() {
        return Module.get_names(this.file_id, this.path, false);
    }
    *values() {
        for (let name of this.keys()) {
            yield this.get(name);
        }
        return;
    }
    *items() {
        for (let name of this.keys()) {
            yield [name, this.get(name)];
        }
        return;
    }
    get_type(obj_path) {
        return Module.get_type(this.file_id, obj_path);
    }
    get_link(obj_path) {
        return Module.get_symbolic_link(this.file_id, obj_path);
    }
    get_external_link(obj_path) {
        return Module.get_external_link(this.file_id, obj_path);
    }
    get(obj_name) {
        let fullpath = (/^\//.test(obj_name)) ? obj_name : this.path + "/" + obj_name;
        fullpath = normalizePath(fullpath);
        let type = this.get_type(fullpath);
        if (type === Module.H5G_GROUP) {
            return new Group(this.file_id, fullpath);
        }
        else if (type === Module.H5G_DATASET) {
            return new Dataset(this.file_id, fullpath);
        }
        else if (type === Module.H5G_LINK) {
            // if get_type succeeds, then get_link must as well
            let target = this.get_link(fullpath);
            return new BrokenSoftLink(target);
        }
        else if (type === Module.H5G_UDLINK) {
            // if get_type succeeds, then get_external_link must as well
            let { filename, obj_path } = this.get_external_link(fullpath);
            return new ExternalLink(filename, obj_path);
        }
        else if (type === Module.H5G_TYPE) {
            return new Datatype(this.file_id, fullpath);
        }
        // unknown type or object not found
        return null;
    }
    create_group(name, track_order = false) {
        Module.create_group(this.file_id, this.path + "/" + name, track_order);
        return this.get(name);
    }
    create_dataset(args) {
        const { name, data, shape, dtype, maxshape, chunks, compression, compression_opts, track_order } = args;
        let metadata = dtype ? dtype_to_metadata(dtype) : guess_metadata(data);
        if (compression && !chunks) {
            throw new Error("cannot specify compression without chunks");
        }
        if (compression && metadata.vlen) {
            throw new Error("cannot specify compression with VLEN data");
        }
        if (!metadata.littleEndian) {
            throw new Error("create_dataset with big-endian dtype is not supported");
        }
        const { data: prepared_data, shape: guessed_shape, maxshape: final_maxshape } = prepare_data(data, metadata, shape, maxshape);
        const final_shape = shape?.map(BigInt) ?? guessed_shape;
        const final_chunks = (chunks) ? chunks.map(BigInt) : null;
        let compression_opts_out;
        let compression_id = 0;
        if (compression && typeof compression === 'number') {
            if (typeof compression_opts === 'undefined') {
                // handle special case where no opts are given, use 'gzip'
                // and numeric value of compression as compression_opts
                compression_id = 1;
                compression_opts_out = [compression];
            }
            else {
                // promote single number to opts list.
                compression_id = compression;
                compression_opts_out = (typeof compression_opts === 'number') ? [compression_opts] : compression_opts;
            }
        }
        else if (compression === 'gzip') {
            compression_id = 1;
            if (compression_opts === undefined) {
                compression_opts_out = [4]; // default compression level
            }
            else {
                compression_opts_out = (typeof compression_opts === 'number') ? [compression_opts] : compression_opts;
            }
        }
        else {
            compression_id = 0;
            compression_opts_out = [];
        }
        const ds_metadata = {
            ...metadata,
            shape: final_shape,
            maxshape: final_maxshape,
            chunks: final_chunks,
            compression: compression_id,
            compression_opts: compression_opts_out,
            track_order: track_order ?? false
        };
        if (metadata.vlen) {
            Module.create_vlen_str_dataset(this.file_id, this.path + "/" + name, prepared_data, ds_metadata);
        }
        else {
            let data_ptr = check_malloc(prepared_data.byteLength);
            try {
                Module.HEAPU8.set(prepared_data, data_ptr);
                Module.create_dataset(this.file_id, this.path + "/" + name, BigInt(data_ptr), ds_metadata);
            }
            finally {
                Module._free(data_ptr);
            }
        }
        return this.get(name);
    }
    create_soft_link(target, name) {
        // create a soft link in this group named *name* to *target* (absolute path)
        const link_name = this.path + '/' + name;
        return Module.create_soft_link(this.file_id, target, link_name);
    }
    create_hard_link(target, name) {
        // create a hard link in this group named *name* to *target* (absolute path)
        const link_name = this.path + '/' + name;
        return Module.create_hard_link(this.file_id, target, link_name);
    }
    create_external_link(file_name, target, name) {
        // create a soft link in this group named *name* to *target* (absolute path)
        const link_name = this.path + '/' + name;
        return Module.create_external_link(this.file_id, file_name, target, link_name);
    }
    toString() {
        return `Group(file_id=${this.file_id}, path=${this.path})`;
    }
    paths() {
        // get all paths below this group in the tree
        return Module.get_names(this.file_id, this.path, true);
    }
}
export class File extends Group {
    constructor(filename, mode = "r", optionsOrTrackOrder = false) {
        let track_order;
        let libver;
        if (optionsOrTrackOrder !== null && typeof optionsOrTrackOrder === "object") {
            track_order = optionsOrTrackOrder.track_order ?? false;
            libver = optionsOrTrackOrder.libver;
        }
        else {
            track_order = optionsOrTrackOrder ?? false;
        }
        const access_mode = ACCESS_MODES[mode];
        const h5_mode = Module[access_mode];
        // Parse libver into numeric bounds
        let libver_low = -1;
        let libver_high = -1;
        if (libver) {
            if (Array.isArray(libver)) {
                // Tuple: [low, high]
                libver_low = parseLibverString(libver[0]);
                libver_high = parseLibverString(libver[1]);
            }
            else {
                // Single string: both bounds set to same value
                const ver = parseLibverString(libver);
                libver_low = ver;
                libver_high = ver;
            }
        }
        const file_id = Module.open(filename, h5_mode, track_order, libver_low, libver_high);
        super(file_id, "/");
        this.filename = filename;
        this.mode = mode;
    }
    get libver() {
        const bounds = Module.get_libver_bounds(this.file_id);
        return [convertToLibverString(bounds.low), convertToLibverString(bounds.high)];
    }
    flush() {
        Module.flush(this.file_id);
    }
    start_swmr_write() {
        return Module.H5Fstart_swmr_write(this.file_id);
    }
    close() {
        return Module.close_file(this.file_id);
    }
}
const calculateHyperslabParams = (shape, ranges) => {
    const strides = shape.map((s, i) => BigInt(ranges?.[i]?.[2] ?? 1));
    const count = shape.map((s, i) => {
        const range_upper = ranges?.[i]?.[1] ?? s;
        const range_lower = ranges?.[i]?.[0] ?? 0;
        const high = (range_upper < s) ? range_upper : s;
        const low = (range_lower > 0) ? range_lower : 0;
        const N = BigInt(high - low);
        const st = strides[i];
        return BigInt(N / st + ((N % st) + st - 1n) / st);
    });
    const offset = shape.map((s, i) => {
        const range_lower = ranges?.[i]?.[0] ?? 0;
        const low = (range_lower > 0) ? range_lower : 0;
        return BigInt((s < low) ? s : low);
    });
    // reurn BigInt arrays, to match inputs of Module functions
    return { strides, count, offset };
};
export class Dataset extends HasAttrs {
    constructor(file_id, path) {
        super();
        this.path = path;
        this.file_id = file_id;
        this.type = OBJECT_TYPE.DATASET;
    }
    refresh() {
        const status = Module.refresh_dataset(this.file_id, this.path);
        if (status < 0) {
            throw new Error(`Could not refresh. Error code: ${status}`);
        }
        delete this._metadata;
    }
    get metadata() {
        if (typeof this._metadata === "undefined") {
            this._metadata = Module.get_dataset_metadata(this.file_id, this.path);
        }
        return this._metadata;
    }
    get dtype() {
        return metadata_to_dtype(this.metadata);
    }
    get shape() {
        return this.metadata.shape;
    }
    get filters() {
        return Module.get_dataset_filters(this.file_id, this.path);
    }
    get value() {
        return this._value_getter(false);
    }
    get json_value() {
        return this._value_getter(true);
    }
    slice(ranges) {
        // interpret ranges as [start, stop], with one per dim.
        const metadata = this.metadata;
        // if auto_refresh is on, getting the metadata has triggered a refresh of the dataset_id;
        const { shape } = metadata;
        if (!shape) {
            return null;
        }
        const { strides, count, offset } = calculateHyperslabParams(shape, ranges);
        const total_size = count.reduce((previous, current) => current * previous, 1n);
        const nbytes = metadata.size * Number(total_size);
        const data_ptr = check_malloc(nbytes);
        let processed;
        try {
            Module.get_dataset_data(this.file_id, this.path, count, offset, strides, BigInt(data_ptr));
            let data = Module.HEAPU8.slice(data_ptr, data_ptr + nbytes);
            processed = process_data(data, metadata, false);
        }
        finally {
            if (metadata.vlen || metadata.type === Module.H5T_class_t.H5T_VLEN.value) {
                Module.reclaim_vlen_memory(this.file_id, this.path, "", BigInt(data_ptr));
            }
            Module._free(data_ptr);
        }
        return processed;
    }
    write_slice(ranges, data) {
        // interpret ranges as [start, stop], with one per dim.
        let metadata = this.metadata;
        if (!metadata.shape) {
            throw new Error("cannot write to a slice of an empty dataset");
        }
        if (metadata.vlen) {
            throw new Error("writing to a slice of vlen dtype is not implemented");
        }
        const { shape } = metadata;
        // if auto_refresh is on, getting the metadata has triggered a refresh of the dataset_id;
        const { strides, count, offset } = calculateHyperslabParams(shape, ranges);
        const { data: prepared_data } = prepare_data(data, metadata, count);
        let data_ptr = check_malloc(prepared_data.byteLength);
        Module.HEAPU8.set(prepared_data, data_ptr);
        try {
            Module.set_dataset_data(this.file_id, this.path, count, offset, strides, BigInt(data_ptr));
        }
        finally {
            Module._free(data_ptr);
        }
    }
    create_region_reference(ranges) {
        const metadata = this.metadata;
        if (!metadata.shape) {
            throw new Error("cannot create region reference from empty dataset");
        }
        // interpret ranges as [start, stop], with one per dim.
        const { shape } = metadata;
        const { strides, count, offset } = calculateHyperslabParams(shape, ranges);
        const ref_data = Module.create_region_reference(this.file_id, this.path, count, offset, strides);
        return new RegionReference(ref_data);
    }
    to_array() {
        const { json_value, metadata } = this;
        const { shape } = metadata;
        if (!isIterable(json_value) || typeof json_value === "string" || shape === null) {
            return json_value;
        }
        let nested = create_nested_array(json_value, shape);
        return nested;
    }
    resize(new_shape) {
        const result = Module.resize_dataset(this.file_id, this.path, new_shape.map(BigInt));
        // reset metadata, pull from file on next read.
        this._metadata = undefined;
        return result;
    }
    make_scale(scale_name = "") {
        // convert dataset to dimension scale
        Module.set_scale(this.file_id, this.path, scale_name);
    }
    attach_scale(index, scale_dset_path) {
        // attach a dimension scale to the "index" dimension of this dataset
        Module.attach_scale(this.file_id, this.path, scale_dset_path, index);
    }
    detach_scale(index, scale_dset_path) {
        // detach a dimension scale from the "index" dimension of this dataset
        Module.detach_scale(this.file_id, this.path, scale_dset_path, index);
    }
    get_attached_scales(index) {
        // get full paths to all datasets that are attached as dimension scales
        // to the specified dimension (at "index") of this dataset.
        return Module.get_attached_scales(this.file_id, this.path, index);
    }
    get_scale_name() {
        // if this dataset is a dimension scale, returns name as string
        // (returns empty string if no name defined, but it is a dimension scale)
        // else returns null if it is not set as a dimension scale
        return Module.get_scale_name(this.file_id, this.path);
    }
    set_dimension_label(index, label) {
        // label dimension at "index" of this dataset with string "label"
        Module.set_dimension_label(this.file_id, this.path, index, label);
    }
    get_dimension_labels() {
        // fetch labels for all dimensions of this dataset (null if label not defined)
        return Module.get_dimension_labels(this.file_id, this.path);
    }
    _value_getter(json_compatible = false) {
        let metadata = this.metadata;
        if (!metadata.shape) {
            return null;
        }
        // if auto_refresh is on, getting the metadata has triggered a refresh of the dataset_id;
        let nbytes = metadata.size * metadata.total_size;
        let data_ptr = check_malloc(nbytes);
        let processed;
        try {
            Module.get_dataset_data(this.file_id, this.path, null, null, null, BigInt(data_ptr));
            let data = Module.HEAPU8.slice(data_ptr, data_ptr + nbytes);
            processed = process_data(data, metadata, json_compatible);
        }
        finally {
            if (metadata.vlen) {
                Module.reclaim_vlen_memory(this.file_id, this.path, "", BigInt(data_ptr));
            }
            Module._free(data_ptr);
        }
        return processed;
    }
}
export class DatasetRegion {
    constructor(source_dataset, region_reference) {
        this.source_dataset = source_dataset;
        this.region_reference = region_reference;
    }
    get metadata() {
        if (typeof this._metadata === "undefined") {
            this._metadata = Module.get_region_metadata(this.source_dataset.file_id, this.region_reference.ref_data);
        }
        return this._metadata;
    }
    get value() {
        return this._value_getter(false);
    }
    _value_getter(json_compatible = false) {
        let metadata = this.metadata;
        if (!metadata.shape) {
            return null;
        }
        // if auto_refresh is on, getting the metadata has triggered a refresh of the dataset_id;
        let nbytes = metadata.size * metadata.total_size;
        let data_ptr = check_malloc(nbytes);
        let processed;
        try {
            Module.get_region_data(this.source_dataset.file_id, this.region_reference.ref_data, BigInt(data_ptr));
            let data = Module.HEAPU8.slice(data_ptr, data_ptr + nbytes);
            processed = process_data(data, metadata, json_compatible);
        }
        finally {
            if (metadata.vlen) {
                Module.reclaim_vlen_memory(this.source_dataset.file_id, this.source_dataset.path, "", BigInt(data_ptr));
            }
            Module._free(data_ptr);
        }
        return processed;
    }
}
function create_nested_array(value, shape) {
    // check that shapes match:
    const total_length = value.length;
    const dims_product = shape.reduce((previous, current) => (previous * current), 1);
    if (total_length !== dims_product) {
        console.warn(`shape product: ${dims_product} does not match length of flattened array: ${total_length}`);
    }
    // Get reshaped output:
    let output = value;
    const subdims = shape.slice(1).reverse();
    for (let dim of subdims) {
        // in each pass, replace input with array of slices of input
        const new_output = [];
        const { length } = output;
        let cursor = 0;
        while (cursor < length) {
            new_output.push(output.slice(cursor, cursor += dim));
        }
        output = new_output;
    }
    return output;
}
export const h5wasm = {
    File,
    Group,
    Dataset,
    Datatype,
    DatasetRegion,
    ready,
    ACCESS_MODES
};
export default h5wasm;
