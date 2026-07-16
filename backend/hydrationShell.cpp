#include "HydrationShell.hpp"

#include <cmath>




std::vector<std::array<double,3>>
HydrationShell::directions =
{
    {{ 1.0, 0.0, 0.0 }},
    {{-1.0, 0.0, 0.0 }},
    {{ 0.0, 1.0, 0.0 }},
    {{ 0.0,-1.0, 0.0 }},
    {{ 0.0, 0.0, 1.0 }},
    {{ 0.0, 0.0,-1.0 }},

    {{ 0.57735027, 0.57735027, 0.57735027 }},
    {{-0.57735027, 0.57735027, 0.57735027 }},
    {{ 0.57735027,-0.57735027, 0.57735027 }},
    {{ 0.57735027, 0.57735027,-0.57735027 }},
    {{-0.57735027,-0.57735027, 0.57735027 }},
    {{-0.57735027, 0.57735027,-0.57735027 }},
    {{ 0.57735027,-0.57735027,-0.57735027 }},
    {{-0.57735027,-0.57735027,-0.57735027 }}
};



std::vector<HydrationPoint>
HydrationShell::generate(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table
) const
{
    std::vector<HydrationPoint> shell;



    for (const auto& atom : atoms)
    {
        if (!table.contains(atom.Z))
            continue;


        const auto& ff =
            table[atom.Z];


        const double shellRadius =
            ff.atomic_radius +
            HydrationShell::WATER_RADIUS;


        for (const auto& dir : directions)
        {
            HydrationPoint p;

            p.x =
                atom.x +
                shellRadius * dir[0];

            p.y =
                atom.y +
                shellRadius * dir[1];

            p.z =
                atom.z +
                shellRadius * dir[2];


            shell.push_back(p);
        }
    }



    return shell;
}