## MODIFIED Requirements

### Requirement: Time-step slider

The interface SHALL provide a horizontal range slider that scrubs through available time steps for the selected file. The slider SHALL display the current time step index. Moving the slider SHALL load the corresponding HDF5 file and update all panels. Slider input SHALL be debounced: after the user stops dragging for 150 ms, a single fetch request SHALL be issued. If a prior fetch is still in flight when a new debounced request fires, the prior request SHALL be aborted via AbortController. The slider label SHALL update immediately on input (before debounce) for visual responsiveness.

#### Scenario: Slider selects time step

- **WHEN** the user drags the slider to position 5
- **THEN** the application SHALL load the file at index 5 from the sorted file list

#### Scenario: Slider debounces rapid scrubbing

- **WHEN** the user rapidly scrubs the slider from step 0 to step 20
- **THEN** only one fetch request SHALL be issued (for step 20) after the 150 ms debounce window
- **AND** the slider label SHALL show the current position continuously during scrubbing

#### Scenario: In-flight request aborted on new selection

- **WHEN** the user selects step 10 and then immediately selects step 15 before step 10 finishes loading
- **THEN** the fetch for step 10 SHALL be aborted
- **AND** only step 15 SHALL be fetched and rendered

### Requirement: File selector control

The interface SHALL provide a dropdown listing all available HDF5 files. Selecting a file SHALL load its parameters and update all panels. If a prior file load is still in flight when a new file is selected, the prior request SHALL be aborted via AbortController. Each controller instance SHALL maintain a single AbortController that is aborted and replaced on each new fetch sequence (file change or slider change).

#### Scenario: File selection triggers load

- **WHEN** the user selects "Output00005.gzip.h5" from the dropdown
- **THEN** the application SHALL fetch params and field data for that file

#### Scenario: File change aborts prior requests

- **WHEN** the user selects file A and then immediately selects file B
- **THEN** all in-flight requests for file A SHALL be aborted
- **AND** only file B data SHALL be loaded

### Requirement: Memory cleanup on panel destroy

When a PanelView is destroyed, the application SHALL:
1. Call `destroy()` on the PanelView's FieldCanvas, removing all model event listeners
2. Call `destroy()` on the PanelView's ColourBar, removing all model event listeners
3. Call `destroy()` on the PanelView's TitleInput, removing all model event listeners
4. Remove the panel's DOM subtree from the document
5. Set the panel's `fieldData` to `null` in the model (or move the full PanelState to stash)

This ensures that detached canvases, large typed arrays, and stale event listeners do not accumulate in memory.

#### Scenario: Destroy removes all listeners

- **WHEN** a panel is destroyed during layout shrink
- **THEN** the model's internal listener count for panel-scoped events SHALL not grow

#### Scenario: Field data freed on non-stashed destroy

- **WHEN** a panel is destroyed and NOT moved to stash
- **THEN** its `fieldData` SHALL be set to `null` and eligible for garbage collection

### Requirement: Initial render timing

All panels SHALL display their field data once the initial data load completes. The application SHALL NOT rely solely on ResizeObserver to trigger the first render. After data arrives for a panel, the application SHALL explicitly trigger a render pass for that panel's canvas and colour bar. If the canvas has zero dimensions at the time data arrives, the render SHALL be deferred until the next ResizeObserver callback reports non-zero dimensions.

#### Scenario: First load renders immediately

- **WHEN** the application starts and loads the first file
- **THEN** all panels SHALL display rendered field data without requiring a manual resize

#### Scenario: Data arrives before layout settled

- **WHEN** field data arrives but the canvas element has zero width
- **THEN** the render SHALL be deferred until the ResizeObserver reports non-zero dimensions
- **AND** the deferred render SHALL proceed automatically without user interaction

### Requirement: Theme toggle (light/dark)

The application SHALL support two visual themes: dark (default) and light. A toggle button in the header SHALL switch between them. All colours used in canvas rendering (background, tick marks, axis labels, frame, null pixels) SHALL be read from CSS custom properties on the document element. Switching themes SHALL update the `data-theme` attribute on `<html>`, fire a `theme-changed` event on the model, and re-render all canvases immediately.

#### Scenario: Toggle from dark to light

- **WHEN** the user clicks the theme toggle
- **THEN** the `data-theme` attribute SHALL be set to `"light"`
- **AND** all canvases SHALL re-render with light-theme colours (white background, dark ticks)

#### Scenario: Canvas reads theme from CSS

- **WHEN** the canvas renders
- **THEN** background colour, tick colour, label colour, frame colour, and null-pixel colour SHALL be read via `getComputedStyle()` from CSS custom properties

### Requirement: Spatial zoom/crop controls

Each panel SHALL provide x-min, x-max, z-min, z-max text inputs for spatially cropping the rendered view. Values SHALL be entered in the current display unit (km or m). The model SHALL store view bounds in SI metres. When view bounds are set, the renderer SHALL compute crop indices from coordinate arrays and render only the cropped sub-region. The canvas aspect ratio SHALL adapt to the cropped dimensions. A Reset button SHALL appear only when zoom is active; clicking it clears the crop. An "Apply to all" button SHALL appear only after a local zoom change; clicking it copies the zoom to all other panels and then hides.

#### Scenario: Crop to sub-region

- **WHEN** the user enters x: -50 to 50 and z: -80 to 0 (in km)
- **THEN** the canvas SHALL render only the sub-region within those spatial bounds
- **AND** axis ticks SHALL reflect the cropped coordinate range

#### Scenario: Reset button visibility

- **WHEN** view bounds are null (full extent)
- **THEN** the Reset button SHALL be hidden
- **WHEN** view bounds are set
- **THEN** the Reset button SHALL be visible

#### Scenario: Apply to all panels

- **WHEN** the user clicks "Apply to all"
- **THEN** all other panels SHALL adopt the same view bounds
- **AND** the "Apply to all" button SHALL hide on all panels

#### Scenario: Spatial unit conversion

- **WHEN** the user switches from km to m
- **THEN** zoom input values SHALL update to reflect the new unit

### Requirement: Dataset (directory) selector

When the server is started with a root directory, it SHALL scan all subdirectories (one level deep) for folders containing `Output*.gzip.h5` files. The application SHALL display a Dataset dropdown above the File selector listing all discovered datasets with their file counts. Selecting a dataset SHALL switch the server's active data directory, reload the file list, and load the first time step. The dropdown SHALL be hidden when only one dataset is found.

#### Scenario: Multiple datasets found

- **WHEN** the root directory contains subdirectories `sim-v1/` (10 files) and `sim-v2/` (20 files)
- **THEN** the Dataset dropdown SHALL list both with file counts
- **AND** the first dataset SHALL be selected by default

#### Scenario: Single dataset hides dropdown

- **WHEN** only one subdirectory contains HDF5 files
- **THEN** the Dataset dropdown SHALL be hidden

#### Scenario: Switch dataset reloads files

- **WHEN** the user selects a different dataset
- **THEN** the file list SHALL reload from the new directory
- **AND** the first time step SHALL be loaded and displayed
