@compute @workgroup_size(64)
fn main(@builtin(global_invocation_id) id : vec3<u32>)
{
    let qIndex = id.x;
    if (qIndex >= uni.qCount) { return; }

    let q = qValues[qIndex];
    var I : f32 = 0.0;

    // Precompute f_i(q) for all atoms ONCE per thread
    // (still O(N), but removes inner-loop overhead)

    for (var i = 0u; i < uni.atomCount; i = i + 1u)
    {
        let fi = formFactors[atoms[i].Z * uni.qCount + qIndex];

        let ri = atoms[i].pos;

        for (var j = 0u; j < uni.atomCount; j = j + 1u)
        {
            let rj = atoms[j].pos;

            // ---- OPTIMISED DISTANCE (NO distance()) ----
            let dx = ri.x - rj.x;
            let dy = ri.y - rj.y;
            let dz = ri.z - rj.z;

            let r2 = dx*dx + dy*dy + dz*dz;
            let r = sqrt(r2);

            let x = q * r;

            let fj = formFactors[atoms[j].Z * uni.qCount + qIndex];

            I = I + fi * fj * sinc(x);
        }
    }

    intensity[qIndex] = I;
}