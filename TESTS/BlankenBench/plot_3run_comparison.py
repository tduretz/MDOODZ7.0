#!/usr/bin/env python3
"""Compare Blankenbach time series: Run 5-8 across resolutions and reseed modes."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import csv

def load_csv(path):
    steps, nu_top, nu_bot, vrms, nb = [], [], [], [], []
    with open(path) as f:
        for line in f:
            if line.startswith('#') or line.startswith('step'):
                continue
            parts = line.strip().split(',')
            if len(parts) < 5:
                continue
            steps.append(int(parts[0]))
            nu_top.append(float(parts[1]))
            nu_bot.append(float(parts[2]))
            vrms.append(float(parts[3]))
            nb.append(int(parts[4]))
    return np.array(steps), np.array(nu_top), np.array(nu_bot), np.array(vrms), np.array(nb)

r5 = load_csv('run5_reference.csv')
r6 = load_csv('run6_reference.csv')
r7 = load_csv('run7_reference.csv')
r8 = load_csv('run8_reference.csv')

Nu_ref = 4.884409
Vrms_ref = 42.864947

fig, axes = plt.subplots(4, 1, figsize=(12, 14), sharex=True)
fig.suptitle('Blankenbach Case 1a — Run 5 / 6 / 7 / 8\n'
             '101² (mode 1,2) · 201² (mode 2) · 401² (mode 2)',
             fontsize=14, fontweight='bold')

colors = {'r5': '#5B9BD5', 'r6': '#ED7D31', 'r7': '#70AD47', 'r8': '#C00000'}
alpha_low = 0.55

# Panel 1: Nu_top
ax = axes[0]
ax.set_title('Nusselt Number (top surface)', fontsize=11)
ax.plot(r5[0], r5[1], color=colors['r5'], alpha=alpha_low, lw=1.0, label='Run 5 (101², mode=1)')
ax.plot(r6[0], r6[1], color=colors['r6'], alpha=alpha_low, lw=1.0, label='Run 6 (101², mode=2)')
ax.plot(r7[0], r7[1], color=colors['r7'], alpha=0.8, lw=1.2, label='Run 7 (201², mode=2)')
ax.plot(r8[0], r8[1], color=colors['r8'], alpha=0.9, lw=1.5, label='Run 8 (401², mode=2)')
ax.axhline(Nu_ref, color='k', ls='--', lw=1.0, label=f'Blankenbach ref = {Nu_ref}')
ax.set_ylabel('Nu_top')
ax.set_ylim(0, 10)
ax.legend(loc='upper right', fontsize=9)

# Panel 2: Nu_bot
ax = axes[1]
ax.set_title('Nusselt Number (bottom surface)', fontsize=11)
ax.plot(r5[0], r5[2], color=colors['r5'], alpha=alpha_low, lw=1.0, label='Run 5')
ax.plot(r6[0], r6[2], color=colors['r6'], alpha=alpha_low, lw=1.0, label='Run 6')
ax.plot(r7[0], r7[2], color=colors['r7'], alpha=0.8, lw=1.2, label='Run 7')
ax.plot(r8[0], r8[2], color=colors['r8'], alpha=0.9, lw=1.5, label='Run 8')
ax.axhline(Nu_ref, color='k', ls='--', lw=1.0, label=f'ref = {Nu_ref}')
ax.set_ylabel('Nu_bot')
ax.set_ylim(0, 6)
ax.legend(loc='lower right', fontsize=9)

# Panel 3: Vrms (clip outliers at 60 for readability)
ax = axes[2]
ax.set_title('RMS Velocity', fontsize=11)
ax.plot(r5[0], np.clip(r5[3], 0, 60), color=colors['r5'], alpha=alpha_low, lw=1.0, label='Run 5')
ax.plot(r6[0], np.clip(r6[3], 0, 60), color=colors['r6'], alpha=alpha_low, lw=1.0, label='Run 6')
ax.plot(r7[0], np.clip(r7[3], 0, 60), color=colors['r7'], alpha=0.8, lw=1.2, label='Run 7')
ax.plot(r8[0], np.clip(r8[3], 0, 60), color=colors['r8'], alpha=0.9, lw=1.5, label='Run 8')
ax.axhline(Vrms_ref, color='k', ls='--', lw=1.0, label=f'ref = {Vrms_ref:.2f}')
ax.set_ylabel('Vrms (non-dim)')
ax.set_ylim(0, 60)
ax.legend(loc='upper right', fontsize=9)

# Panel 4: Particle count
ax = axes[3]
ax.set_title('Particle Count', fontsize=11)
ax.plot(r5[0], r5[4], color=colors['r5'], alpha=alpha_low, lw=1.0, label='Run 5 (101², mode=1)')
ax.plot(r6[0], r6[4], color=colors['r6'], alpha=alpha_low, lw=1.0, label='Run 6 (101², mode=2)')
ax.plot(r7[0], r7[4], color=colors['r7'], alpha=0.8, lw=1.2, label='Run 7 (201², mode=2)')
ax.plot(r8[0], r8[4], color=colors['r8'], alpha=0.9, lw=1.5, label='Run 8 (401², mode=2)')
ax.set_ylabel('Nb_part')
ax.set_xlabel('Step')
ax.legend(loc='right', fontsize=9)

plt.tight_layout()
plt.savefig('blankenbach_run5_vs_run6_vs_run7_vs_run8.png', dpi=150)
print('Saved: blankenbach_run5_vs_run6_vs_run7_vs_run8.png')
