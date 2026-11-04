struct Params {
    atomCount : u32,
    qCount : u32,
    qOffset : u32,
    iTileCount : u32,
};

@group(0) @binding(0)
var<storage, read> positions : array<f32>;

@group(0) @binding(1)
var<storage, read> amplitudes : array<f32>;

@group(0) @binding(2)
var<storage, read_write> output : array<f32>;

@group(0) @binding(3)
var<storage, read> qValues : array<f32>;

@group(0) @binding(4)
var<uniform> params : Params;


// Shared memory for the j dimension
var<workgroup> tileX : array<f32, 64>;
var<workgroup> tileY : array<f32, 64>;
var<workgroup> tileZ : array<f32, 64>;
var<workgroup> tileF : array<f32, 64>;


// Shared memory for reduction over i
var<workgroup> partial : array<f32, 64>;


fn sinc(x : f32) -> f32 {

    if (abs(x) < 1e-6) {
        return 1.0;
    }

    return sin(x) / x;
}


@compute
@workgroup_size(64)
fn main(
    @builtin(workgroup_id)
    workgroupId : vec3<u32>,

    @builtin(local_invocation_id)
    localId : vec3<u32>
)
{
    let qIndex = workgroupId.x;
    let iTile = workgroupId.y;
    let localIndex = localId.x;

    if (qIndex >= params.qCount || iTile >= params.iTileCount) {
        return;
    }

    let q = qValues[qIndex];

    // Amplitudes for this q are stored contiguously.
    let amplitudeOffset = qIndex * params.atomCount;

    // This workgroup handles one block of 64 i atoms.
    let i = iTile * 64u + localIndex;

    var xi : f32 = 0.0;
    var yi : f32 = 0.0;
    var zi : f32 = 0.0;
    var fi : f32 = 0.0;

    if (i < params.atomCount) {
        let ix = i * 3u;

        xi = positions[ix];
        yi = positions[ix + 1u];
        zi = positions[ix + 2u];

        fi = amplitudes[amplitudeOffset + i];
    }

    var total : f32 = 0.0;

    let jTileCount = (params.atomCount + 63u) / 64u;

    // Loop over all j tiles.
    for (
        var jTile : u32 = 0u;
        jTile < jTileCount;
        jTile++
    ) {
        let j = jTile * 64u + localIndex;

        // Load j atoms into shared memory.
        if (j < params.atomCount) {
            let jx = j * 3u;

            tileX[localIndex] = positions[jx];
            tileY[localIndex] = positions[jx + 1u];
            tileZ[localIndex] = positions[jx + 2u];

            tileF[localIndex] = amplitudes[amplitudeOffset + j];
        } else {
            tileX[localIndex] = 0.0;
            tileY[localIndex] = 0.0;
            tileZ[localIndex] = 0.0;
            tileF[localIndex] = 0.0;
        }

        workgroupBarrier();

        // Calculate this i atom against all j atoms in the tile.
        for (var k : u32 = 0u; k < 64u; k++) {
            let jIndex = jTile * 64u + k;

            if (i < params.atomCount && jIndex < params.atomCount) {
                let dx = xi - tileX[k];
                let dy = yi - tileY[k];
                let dz = zi - tileZ[k];

                let r = sqrt(
                    dx * dx +
                    dy * dy +
                    dz * dz
                );

                let qr = q * r;

                total += fi * tileF[k] * sinc(qr);
            }
        }

        workgroupBarrier();
    }

    // Reduce the 64 i contributions within this block.
    partial[localIndex] = total;

    workgroupBarrier();

    var stride : u32 = 32u;

    while (stride > 0u) {
        if (localIndex < stride) {
            partial[localIndex] += partial[localIndex + stride];
        }

        workgroupBarrier();

        stride = stride / 2u;
    }

    // Each workgroup writes its own partial result.
    if (localIndex == 0u) {
        let outputIndex = qIndex * params.iTileCount + iTile;
        output[outputIndex] = partial[0];
    }
}