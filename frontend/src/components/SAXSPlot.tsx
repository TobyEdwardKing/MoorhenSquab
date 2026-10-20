import {
    Chart as ChartJS,
    CategoryScale,
    LinearScale,
    PointElement,
    LineElement,
    Title,
    Tooltip,
    Legend,
} from "chart.js";

import { Scatter } from "react-chartjs-2";

ChartJS.register(
    CategoryScale,
    LinearScale,
    PointElement,
    LineElement,
    Title,
    Tooltip,
    Legend
);

interface SAXSPlotProps {
    experimentalQ: number[];
    experimentalIntensity: number[];

    calculatedQ: number[];
    calculatedIntensity: number[];

    scaledIntensity?: number[];
}

export function SAXSPlot({
    experimentalQ,
    experimentalIntensity,
    calculatedQ,
    calculatedIntensity,
    scaledIntensity,
}: SAXSPlotProps) {

const datasets = [
    {
        label: "Experimental",
        data: experimentalQ.map((q, i) => ({
            x: q,
            y: experimentalIntensity[i],
        })),
        showLine: false,
        pointRadius: 2,
        backgroundColor: "blue",
        borderColor: "blue",
    },
];

if (scaledIntensity) {
    datasets.push({
        label: "Calculated",
        data: calculatedQ.map((q, i) => ({
            x: q,
            y: scaledIntensity[i],
        })),
        showLine: true,
        pointRadius: 0,
        backgroundColor: "orange",
        borderColor: "orange",
    });
}

    const data = {
        datasets,
    };

    const options = {
        responsive: true,

        plugins: {
            legend: {
                display: true,
            },

            title: {
                display: true,
                text: "SAXS profile",
            },
        },

        scales: {
            x: {
                type: "linear" as const,

                title: {
                    display: true,
                    text: "q (Å⁻¹)",
                },
            },

            y: {
                title: {
                    display: true,
                    text: "Intensity",
                },
            },
        },
    };

    return (
        <div style={{ width: "100%", maxWidth: "1000px" }}>
            <Scatter
                data={data}
                options={options}
            />
        </div>
    );
}