#include "hydrationShell.hpp"

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


static bool point_is_buried(
    double x,
    double y,
    double z,
    const AtomData& sourceAtom,
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table)
{
    for (const auto& other : atoms)
    {
        // Don't reject against the atom that generated this point.
        if (&other == &sourceAtom)
            continue;

        if (!table.contains_element(other.Z))
            continue;

        const auto& otherFF =
            table.get_element(other.Z);

        const double dx = x - other.x;
        const double dy = y - other.y;
        const double dz = z - other.z;

        const double distanceSquared =
            dx*dx + dy*dy + dz*dz;

        const double radius =
            otherFF.atomic_radius;

        if (distanceSquared < radius * radius)
            return true;
    }

    return false;
}

std::vector<HydrationPoint>
HydrationShell::generate(
    const std::vector<AtomData>& atoms,
    const FormFactorTable& table
) const
{
    std::vector<HydrationPoint> shell;



    for (const auto& atom : atoms)
    {
        if (!table.contains_element(atom.Z))
            continue;


        const auto& ff =
            table.get_element(atom.Z);


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

            if (point_is_buried(
                    p.x,
                    p.y,
                    p.z,
                    atom,
                    atoms,
                    table))
            {
                continue;
            }

            shell.push_back(p);
        }
    }



    return shell;
}