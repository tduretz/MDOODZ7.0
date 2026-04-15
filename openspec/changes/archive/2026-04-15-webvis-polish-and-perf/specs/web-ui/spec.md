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
