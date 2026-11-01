"""Build the exact conjugate-position comparison. Python 3.12 + Matplotlib.

Usage: python3.12 build_conjugate_diagrams.py [--preview /absolute/path.png]
The SVG is a schematic of variable roles, not a metric on dimensional axes.
"""
import argparse
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Circle


def main():
    if sys.version_info[:2] != (3, 12):
        raise RuntimeError('Use Python 3.12')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--preview', type=Path)
    args = parser.parse_args()
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'svg.fonttype': 'path',
                         'svg.hashsalt': 'nsr2-conjugate-diagrams-v1'})
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.6))
    fig.patch.set_facecolor('#f8fafc')
    rows = [
        ('THERMODYNAMICS', [r'$T$', r'$V$', r'$S$', r'$P$'],
         r'$U(S,V)$', r'$dU=T\,dS-P\,dV$'),
        ('ELECTROMAGNETISM', [r'$E$', r'$B$', r'$D$', r'$-H$'],
         r'$u(D,B)$', r'$du=E\,dD+H\,dB$'),
        ('CHARGE–FLUX MODEL', [r'$v$', r'$\Phi$', r'$Q$', r'$-I$'],
         r'$U(Q,\Phi)$', r'$dU=v\,dQ+I\,d\Phi$'),
    ]
    for ax, (title, labels, potential, law) in zip(axes, rows):
        ax.set_aspect('equal')
        ax.set_xlim(-1.7, 1.7)
        ax.set_ylim(-1.85, 1.85)
        ax.axis('off')
        ax.add_patch(Circle((0, 0), 1, fill=False, edgecolor='#94a3b8', linewidth=1.5))
        ax.plot([0, 0], [-0.91, 0.91], color='#0f766e', linewidth=2)
        ax.plot([-0.91, 0.91], [0, 0], color='#b45309', linewidth=2)
        for (x, y), label, color in zip(
                [(0, 1.18), (1.25, 0), (0, -1.18), (-1.25, 0)],
                labels, ['#0f766e', '#b45309', '#0f766e', '#b45309']):
            ax.text(x, y, label, ha='center', va='center', fontsize=23, color=color)
        ax.text(0, 0, potential, ha='center', va='center', fontsize=14,
                bbox={'facecolor': '#f8fafc', 'edgecolor': 'none', 'pad': 5})
        ax.text(0, 1.64, title, ha='center', fontsize=10.5, fontweight='bold', color='#0f172a')
        ax.text(0, -1.64, law, ha='center', fontsize=13, color='#334155')
    fig.suptitle('The same conjugate-pair pattern carries the lambda / z definitions',
                 fontsize=16, fontweight='bold', y=0.98, color='#0f172a')
    fig.text(0.5, 0.035,
             'Schematic positions. Two independent state variables per diagram; signs follow the work differential.',
             ha='center', fontsize=10, color='#475569')
    fig.subplots_adjust(left=0.02, right=0.98, top=0.89, bottom=0.11, wspace=0.08)
    target = Path(__file__).with_name('conjugate_diagrams.svg')
    fig.savefig(target,
                metadata={'Date': None, 'Creator': 'NSR2 Python 3.12 figure builder'})
    target.write_text('\n'.join(line.rstrip() for line in target.read_text().splitlines())+'\n')
    if args.preview:
        fig.savefig(args.preview, dpi=160, facecolor=fig.get_facecolor())
    plt.close(fig)


if __name__ == '__main__':
    main()
