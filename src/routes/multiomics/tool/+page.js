// The analysis workbench relies entirely on browser file APIs and intensive
// interactive Svelte components; its static SSR markup is not part of the
// portable or scientific contract. Kit 3 hydration of this large form triggers
// an HTML hierarchy mismatch in production preview. Keep the workbench a
// client-rendered route while preserving SSR/prerender for lessons and shell.
export const ssr = false;
