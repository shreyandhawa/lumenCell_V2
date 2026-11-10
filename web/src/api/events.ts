import { useEffect, useState } from "react";
import { jobEventsUrl } from "@/api/client";

/** One server-sent progress event (mirrors the server's JobEvent). */
export interface JobEvent {
  job_id: string;
  state: "queued" | "running" | "done" | "failed" | "cancelled";
  stage: string | null;
  i: number;
  n: number;
  fraction: number;
  message: string | null;
  elapsed: number;
  eta: number | null;
}

const TERMINAL = new Set(["done", "failed", "cancelled"]);

/** Follow a job's progress over SSE until it finishes. */
export function useJobEvents(jobId: string | null | undefined): JobEvent | null {
  const [event, setEvent] = useState<JobEvent | null>(null);
  useEffect(() => {
    // EventSource is missing in some environments (jsdom); callers fall back to polled status.
    if (!jobId || typeof EventSource === "undefined") return;
    const source = new EventSource(jobEventsUrl(jobId));
    source.onmessage = (msg: MessageEvent<string>) => {
      const next = JSON.parse(msg.data) as JobEvent;
      setEvent(next);
      if (TERMINAL.has(next.state)) source.close();
    };
    source.onerror = () => source.close();
    return () => {
      source.close();
      setEvent(null);
    };
  }, [jobId]);
  return event;
}
