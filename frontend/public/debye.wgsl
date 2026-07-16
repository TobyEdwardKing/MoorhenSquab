struct Params {
    atomCount : u32,
    q : f32,
};

@group(0) @binding(0)
var<storage, read> positions : array<f32>;

@group(0) @binding(1)
var<storage, read> formFactors : array<f32>;

@group(0) @binding(2)
var<storage, read_write> output : array<f32>;

@group(0) @binding(3)
var<uniform> params : Params;

fn sinc(x : f32) -> f32 {

    if (abs(x) < 1e-6) {
        return 1.0;
    }

    return sin(x) / x;
}

@compute
@workgroup_size(64)
fn main(
    @builtin(global_invocation_id)
    id : vec3<u32>
)
{
    let i = id.x;

    if (i >= params.atomCount) {
        return;
    }

    let ix = i * 3u;

    let xi = positions[ix];
    let yi = positions[ix + 1u];
    let zi = positions[ix + 2u];

    let fi = formFactors[i];

    var total : f32 = 0.0;

    for (var j : u32 = 0u; j < params.atomCount; j++) {

        let jx = j * 3u;

        let dx = xi - positions[jx];
        let dy = yi - positions[jx + 1u];
        let dz = zi - positions[jx + 2u];

        let r = sqrt(dx*dx + dy*dy + dz*dz);


        let qr = params.q * r;

        total +=
            fi *
            formFactors[j] *
            sinc(qr);
        }
    output[i] = total;
}