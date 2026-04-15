// ── TitleInput View ───────────────────────────────────────────────────
// Editable title input with ${...} token dropdown.

const TOKEN_CATALOGUE = [
  { group: 'Metadata', items: [
    { token: 'step',       label: 'Step number' },
    { token: 'time',       label: 'Model time' },
    { token: 'Nx',         label: 'Grid columns' },
    { token: 'Nz',         label: 'Grid rows' },
    { token: 'resolution', label: 'Grid resolution' },
  ]},
  { group: 'Field stats', items: [
    { token: 'min',    label: 'Min' },
    { token: 'max',    label: 'Max' },
    { token: 'mean',   label: 'Mean' },
    { token: 'median', label: 'Median' },
    { token: 'p2',     label: 'Robust min (p2)' },
    { token: 'p98',    label: 'Robust max (p98)' },
    { token: 'field',  label: 'Field name' },
    { token: 'unit',   label: 'Unit' },
  ]},
];

export class TitleInput {
  constructor(model, containerEl, panelState, dispatchTarget) {
    this.model = model;
    this.panelState = panelState;
    this.dispatchTarget = dispatchTarget;

    // Build DOM
    this.row = document.createElement('div');
    this.row.className = 'panel-title-row';

    this.input = document.createElement('input');
    this.input.type = 'text';
    this.input.value = panelState.title;

    this.insertBtn = document.createElement('button');
    this.insertBtn.className = 'insert-btn';
    this.insertBtn.textContent = '⊕';
    this.insertBtn.title = 'Insert data token';

    this.row.appendChild(this.input);
    this.row.appendChild(this.insertBtn);

    containerEl.appendChild(this.row);

    // Dropdown (created lazily)
    this.dropdown = null;

    // Events
    this.input.addEventListener('input', () => this._onInputChange());
    this.input.addEventListener('keydown', e => {
      if (e.key === '$') {
        // Show dropdown after the $ is typed
        setTimeout(() => this._showDropdown(), 0);
      }
      if (e.key === 'Escape' && this.dropdown) {
        this._hideDropdown();
      }
    });
    this.insertBtn.addEventListener('click', () => this._showDropdown());

    // Close dropdown on outside click
    this._outsideClickHandler = (e) => {
      if (this.dropdown && !this.row.contains(e.target)) {
        this._hideDropdown();
      }
    };
    document.addEventListener('click', this._outsideClickHandler, true);
  }

  _onInputChange() {
    this.panelState.title = this.input.value;
    this.dispatchTarget.dispatchEvent(new CustomEvent('ctrl:panel:title', {
      bubbles: true,
      detail: { panelId: this.panelState.id, template: this.input.value }
    }));
  }

  _showDropdown() {
    if (this.dropdown) this._hideDropdown();

    this.dropdown = document.createElement('div');
    this.dropdown.className = 'token-dropdown';

    for (const group of TOKEN_CATALOGUE) {
      const header = document.createElement('div');
      header.className = 'group-header';
      header.textContent = group.group;
      this.dropdown.appendChild(header);

      for (const item of group.items) {
        const el = document.createElement('div');
        el.className = 'token-item';
        el.innerHTML = `<span>${item.label}</span><span class="token-key">\${${item.token}}</span>`;
        el.addEventListener('click', (e) => {
          e.stopPropagation();
          this._insertToken(item.token);
          this._hideDropdown();
        });
        this.dropdown.appendChild(el);
      }
    }

    this.row.appendChild(this.dropdown);
  }

  _hideDropdown() {
    if (this.dropdown) {
      this.dropdown.remove();
      this.dropdown = null;
    }
  }

  _insertToken(tokenName) {
    const text = `\${${tokenName}}`;
    const pos = this.input.selectionStart ?? this.input.value.length;
    const before = this.input.value.slice(0, pos);
    const after = this.input.value.slice(pos);
    this.input.value = before + text + after;
    this.input.selectionStart = this.input.selectionEnd = pos + text.length;
    this.input.focus();
    this._onInputChange();
  }

  destroy() {
    document.removeEventListener('click', this._outsideClickHandler, true);
    this._hideDropdown();
  }
}
