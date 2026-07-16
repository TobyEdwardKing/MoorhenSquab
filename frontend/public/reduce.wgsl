struct Params {
    atomCount : u32,
};

@group(0) @binding(0)
var<storage, read> input : array<f32>;

@group(0) @binding(1)
var<storage, read_write> output : array<f32>;

@group(0) @binding(2)
var<uniform> params : Params;

@compute
@workgroup_size(1)
fn main()
{
    var total : f32 = 0.0;

    for (var i : u32 = 0u; i < params.atomCount; i++)
    {
        total += input[i];
    }

    output[0] = total;
}