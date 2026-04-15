## Why

The WebVis multi-panel feature (shipped in `webvis-multi-panel-titles`) works but has six usability and performance issues reported during testing:

1. **Canvas blank on initial load** — after selecting a field, FieldCanvas sometimes doesn't render because the ResizeObserver fires before field data arrives.
2. **2×1 / 1×2 layouts distort on window resize** — the aspect-ratio logic fails for tall-narrow or wide-short grid cells because `availW` and `availH` compute incorrectly.
3. **Layout round-trip resets panels** — going 2×2 → 1×1 → 2×2 destroys panels 2–4 and re-creates them with blank state instead of preserving them.
4. **Title lives in DOM, not canvas** — the rendered title is a `<div>` overlay; the user wants it baked into the canvas so it appears in screenshot/download.
5. **Slider scrubbing floods requests** — dragging the time slider fires a `selectFile()` for every slider position, creating dozens of parallel fetches that saturate the browser. Needs debounce + AbortController.
6. **Memory leak (622 MB+)** — field data arrays accumulate in PanelState without being freed when new data arrives, and detached event listeners from destroyed PanelViews are never cleaned up.

## What Changes

- **Fix initial render timing**: ensure `FieldCanvas.render()` is called after data arrives, not only on resize
- **Fix 2×1 / 1×2 aspect-ratio math**: use the panel element's actual bounding rect, not a hardcoded `availW - 88`
- **Preserve panel states on layout shrink/grow**: keep removed panels in a stash and restore them when the layout grows back
- **Render title text into the canvas**: draw the interpolated title as a top overlay line on the FieldCanvas, remove the separate DOM `<span>`
- **Debounce slider + abort in-flight requests**: wrap `selectFile()` with a debounce (150 ms) and an `AbortController` that cancels superseded fetches
- **Fix memory leaks**: null out previous `fieldData` before assigning new data, explicitly remove all event listeners in `PanelView.destroy()`, and clean up FieldCanvas / ColourBar / TitleInput references

## Capabilities

### New Capabilities

(none)

### Modified Capabilities
- `field-renderer`: title text now rendered inside the canvas; canvas sizing logic changed
- `multi-panel-layout`: panel stash for layout round-trips; responsive sizing fix
- `web-ui`: slider debounce + fetch abort; memory cleanup

## Impact

- **Files**: `app.mjs`, `controller.mjs`, `model.mjs`, `FieldCanvas.mjs`, `ColourBar.mjs`, `PanelView.mjs`, `TitleInput.mjs`, `style.css`
- **APIs**: no server changes
- **Dependencies**: none new
- **Risk**: low — all changes are client-side UX/perf fixes with no data-format changes
