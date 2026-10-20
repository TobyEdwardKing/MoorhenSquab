import React, { useState } from "react";
import { parseSAXSProfile } from "../utils/SAXSProfileParser";
import { SAXSProfile } from "../types/SAXSProfile";

export interface SAXSProfileLoaderProps {
    onProfileLoaded: (profile: SAXSProfile, fileName: string) => void;
}

export function SAXSProfileLoader({
    onProfileLoaded
}: SAXSProfileLoaderProps) {

    const [loading, setLoading] = useState(false);
    const [error, setError] = useState<string | null>(null);
    const [fileName, setFileName] = useState<string | null>(null);

    const handleFileChange = async (
        event: React.ChangeEvent<HTMLInputElement>
    ) => {

        const file = event.target.files?.[0];

        if (!file) {
            return;
        }

        setLoading(true);
        setError(null);
        setFileName(null);

        try {

            const profile = await parseSAXSProfile(file);

            setFileName(file.name);

            onProfileLoaded(
                profile,
                file.name
            );

        } catch (err) {

            console.error(
                "Failed to load SAXS profile:",
                err
            );

            setError(
                err instanceof Error
                    ? err.message
                    : "Failed to load SAXS profile."
            );

        } finally {

            setLoading(false);
        }
    };

    return (
        <div>
            <label>
                Load SAXS profile
            </label>

            <input
                type="file"
                accept=".dat,.txt,.csv"
                onChange={handleFileChange}
                disabled={loading}
            />

            {loading && (
                <p>
                    Loading SAXS profile...
                </p>
            )}

            {fileName && !loading && (
                <p>
                    Loaded: {fileName}
                </p>
            )}

            {error && (
                <p>
                    Error: {error}
                </p>
            )}
        </div>
    );
}
