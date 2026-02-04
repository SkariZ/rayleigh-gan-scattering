import sys
import numpy as np

from PyQt5 import QtWidgets
from PyQt5.QtCore import Qt

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

from rayleigh_gans.form_factors import (
    q_from_theta,
    F_sphere_q,
    F_thin_spherical_shell_q,
    F_thick_spherical_shell_q,
    F_disk_q,
    P_from_F,
    P_gaussian_coil_q,
    P_guinier_q,
)

# ----------------------------
# UI helpers
# ----------------------------

def make_slider(min_v, max_v, init_v, step=1):
    s = QtWidgets.QSlider(Qt.Horizontal)
    s.setMinimum(int(min_v))
    s.setMaximum(int(max_v))
    s.setSingleStep(int(step))
    s.setValue(int(init_v))
    s.setTracking(True)
    return s


class MplCanvas(FigureCanvas):
    def __init__(self, parent=None):
        fig = Figure(figsize=(6, 4), dpi=140)
        self.ax = fig.add_subplot(111)
        super().__init__(fig)
        self.setParent(parent)


# ----------------------------
# Geometry registry
# ----------------------------

GEOMS = [
    "Sphere",
    "Thin spherical shell",
    "Thick spherical shell",
    "Disk",
    "Gaussian coil (Debye)",
    "Guinier",
]

GEOM_HELP = {
    "Sphere": (
        "Sphere \n"
        "F(q) = 3[sin(qR) - qR cos(qR)]/(qR)^3\n"
        "Parameters: R (radius)\n"
        "Notes: normalized so F(0)=1"
    ),
    "Thin spherical shell": (
        "Thin spherical shell / surface \n"
        "F(q) = sin(qR)/(qR)\n"
        "Parameters: R (shell radius)\n"
        "Notes: this is a 'skin' (infinitesimally thin)"
    ),
    "Thick spherical shell": (
        "Thick spherical shell \n"
        "Computed as difference of two spheres, normalized by shell volume:\n"
        "F_shell(q) = [Ro^3 F_sphere(q,Ro) - Ri^3 F_sphere(q,Ri)]/(Ro^3 - Ri^3)\n"
        "Parameters: R (outer radius), ΔR (thickness) with Ri = R-ΔR\n"
    ),
    "Disk": (
        "Uniform 2D disk\n"
        "F(q_perp) = 2 J1(q_perp R)/(q_perp R)\n"
        "In this viewer we use the quick-look approximation q_perp ≈ |q|.\n"
        "Parameters: R (disk radius)\n"
        "Notes: orientation averaging can be added later."
    ),
    "Gaussian coil (Debye)": (
        "Gaussian polymer coil (Debye intensity)\n"
        "P(q) = 2[exp(-u) + u - 1]/u^2,  u=(q Rg)^2\n"
        "Parameters: Rg (radius of gyration)\n"
        "Notes: already normalized so P(0)=1"
    ),
    "Guinier": (
        "Guinier approximation (small-q intensity)\n"
        "P(q) ≈ exp(-q^2 Rg^2 / 3)\n"
        "Parameters: Rg\n"
        "Validity: typically q*Rg ≲ 1"
    ),
}


def normalize_curve(y: np.ndarray) -> np.ndarray:
    """Normalize for visualization: y / max(|y|)."""
    y = np.asarray(y, dtype=float)
    m = np.nanmax(np.abs(y))
    if not np.isfinite(m) or m <= 0:
        return y
    return y / m


# ----------------------------
# Main window
# ----------------------------

class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Rayleigh–Gans Scattering Viewer (PyQt5)")
        self.resize(1150, 680)

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        layout = QtWidgets.QHBoxLayout(central)

        controls = QtWidgets.QWidget()
        controls.setFixedWidth(420)
        c_layout = QtWidgets.QVBoxLayout(controls)
        c_layout.setAlignment(Qt.AlignTop)

        self.canvas = MplCanvas()
        layout.addWidget(controls)
        layout.addWidget(self.canvas, stretch=1)

        # ---- Title
        title = QtWidgets.QLabel("Controls")
        title.setStyleSheet("font-size: 16px; font-weight: 600;")
        c_layout.addWidget(title)

        # ---- Geometry checklist (multi-select)
        geom_box = QtWidgets.QGroupBox("Geometries (check multiple)")
        geom_layout = QtWidgets.QVBoxLayout(geom_box)

        self.geom_list = QtWidgets.QListWidget()
        for g in GEOMS:
            item = QtWidgets.QListWidgetItem(g)
            item.setFlags(item.flags() | Qt.ItemIsUserCheckable)
            item.setCheckState(Qt.Unchecked)
            self.geom_list.addItem(item)

        # default: Sphere checked
        self.geom_list.item(0).setCheckState(Qt.Checked)

        geom_layout.addWidget(self.geom_list)

        btn_row = QtWidgets.QHBoxLayout()
        self.btn_all = QtWidgets.QPushButton("All")
        self.btn_none = QtWidgets.QPushButton("Clear")
        btn_row.addWidget(self.btn_all)
        btn_row.addWidget(self.btn_none)
        geom_layout.addLayout(btn_row)

        c_layout.addWidget(geom_box)

        # ---- X-axis mode
        self.xmode = QtWidgets.QComboBox()
        self.xmode.addItems(["theta (radians)", "q"])
        c_layout.addWidget(self._row("X-axis", self.xmode))

        # ---- Y mode (including normalization)
        self.ytype = QtWidgets.QComboBox()
        self.ytype.addItems([
            "Re[F(q)]",
            "|F(q)|",
            "|F(q)|²",
        ])
        c_layout.addWidget(self._row("Y", self.ytype))

        # ---- theta range (radians), scaled x100
        self.theta_min = make_slider(0, 150, 1)    # 0..1.50 rad
        self.theta_max = make_slider(1, 314, 314)  # up to pi
        self.theta_n = make_slider(200, 1000, 800, step=10)

        c_layout.addWidget(self._slider_row("theta min [rad]", self.theta_min, scale=100))
        c_layout.addWidget(self._slider_row("theta max [rad]", self.theta_max, scale=100))
        c_layout.addWidget(self._slider_row("N points", self.theta_n, scale=1))

        # ---- q range, scaled x1000
        self.q_min = make_slider(0, 2000, 0)       # 0..2.0
        self.q_max = make_slider(1, 2000, 800)     # 0.8
        c_layout.addWidget(self._slider_row("q min", self.q_min, scale=1000))
        c_layout.addWidget(self._slider_row("q max", self.q_max, scale=1000))

        line = QtWidgets.QFrame()
        line.setFrameShape(QtWidgets.QFrame.HLine)
        c_layout.addWidget(line)

        # ---- Physical parameters
        self.radius = make_slider(1, 1000, 200)          # 1..1000
        self.rg = make_slider(1, 3000, 300)              # keep wide
        self.shell_thickness = make_slider(1, 1000, 50)

        self.wavelength = make_slider(300, 800, 532)     # 300..800
        self.nm = make_slider(100, 200, 133)             # n scaled x100

        c_layout.addWidget(self._slider_row("radius R", self.radius, scale=1))
        c_layout.addWidget(self._slider_row("Rg", self.rg, scale=1))
        c_layout.addWidget(self._slider_row("shell thickness ΔR", self.shell_thickness, scale=1))
        c_layout.addWidget(self._slider_row("wavelength λ", self.wavelength, scale=1))
        c_layout.addWidget(self._slider_row("n_medium", self.nm, scale=100))

        # ---- Log options (log y OFF by default as requested)
        self.logy = QtWidgets.QCheckBox("log y")
        self.logy.setChecked(False)
        c_layout.addWidget(self.logy)

        self.logx = QtWidgets.QCheckBox("log x (only for q)")
        self.logx.setChecked(False)
        c_layout.addWidget(self.logx)

        # ---- Help panel
        help_box = QtWidgets.QGroupBox("Explanation / help")
        help_layout = QtWidgets.QVBoxLayout(help_box)
        self.help_text = QtWidgets.QTextEdit()
        self.help_text.setReadOnly(True)
        self.help_text.setFixedHeight(150)
        help_layout.addWidget(self.help_text)
        c_layout.addWidget(help_box)

        # ---- Update button
        self.btn_update = QtWidgets.QPushButton("Update")
        c_layout.addWidget(self.btn_update)

        # ---- Save button
        self.btn_save = QtWidgets.QPushButton("Save plot (PNG)")
        c_layout.addWidget(self.btn_save)
        self.btn_save.clicked.connect(self.save_plot)

        # Signals
        self.btn_update.clicked.connect(self.update_plot)
        self.btn_all.clicked.connect(self.check_all_geoms)
        self.btn_none.clicked.connect(self.clear_all_geoms)

        self.geom_list.itemChanged.connect(self.on_geom_changed)

        for w in [
            self.xmode, self.ytype,
            self.theta_min, self.theta_max, self.theta_n,
            self.q_min, self.q_max,
            self.radius, self.rg, self.shell_thickness,
            self.wavelength, self.nm,
            self.logx, self.logy,
        ]:
            if isinstance(w, QtWidgets.QSlider):
                w.valueChanged.connect(self.update_plot)
            elif isinstance(w, QtWidgets.QComboBox):
                w.currentIndexChanged.connect(self.update_plot)
            elif isinstance(w, QtWidgets.QCheckBox):
                w.stateChanged.connect(self.update_plot)

        # init help
        self._update_help_text()

        # initial draw
        self.update_plot()

    # ---- UI builders

    def _row(self, label, widget):
        row = QtWidgets.QWidget()
        l = QtWidgets.QHBoxLayout(row)
        l.setContentsMargins(0, 0, 0, 0)
        lab = QtWidgets.QLabel(label)
        lab.setFixedWidth(130)
        l.addWidget(lab)
        l.addWidget(widget, stretch=1)
        return row

    def _slider_row(self, label, slider, scale=1):
        row = QtWidgets.QWidget()
        l = QtWidgets.QVBoxLayout(row)
        l.setContentsMargins(0, 0, 0, 0)

        top = QtWidgets.QWidget()
        ht = QtWidgets.QHBoxLayout(top)
        ht.setContentsMargins(0, 0, 0, 0)

        lab = QtWidgets.QLabel(label)
        val = QtWidgets.QLabel("")
        val.setAlignment(Qt.AlignRight)
        ht.addWidget(lab)
        ht.addWidget(val, stretch=1)

        def sync_label():
            v = slider.value() / scale
            val.setText(f"{v:.4g}" if scale != 1 else f"{int(slider.value())}")

        slider.valueChanged.connect(sync_label)
        sync_label()

        l.addWidget(top)
        l.addWidget(slider)
        return row

    # ---- Geometry selection

    def selected_geometries(self):
        geoms = []
        for i in range(self.geom_list.count()):
            item = self.geom_list.item(i)
            if item.checkState() == Qt.Checked:
                geoms.append(item.text())
        return geoms

    def check_all_geoms(self):
        self.geom_list.blockSignals(True)
        for i in range(self.geom_list.count()):
            self.geom_list.item(i).setCheckState(Qt.Checked)
        self.geom_list.blockSignals(False)
        self._update_help_text()
        self.update_plot()

    def clear_all_geoms(self):
        self.geom_list.blockSignals(True)
        for i in range(self.geom_list.count()):
            self.geom_list.item(i).setCheckState(Qt.Unchecked)
        self.geom_list.blockSignals(False)
        self._update_help_text()
        self.update_plot()

    def on_geom_changed(self, _item):
        self._update_help_text()
        self.update_plot()

    def _update_help_text(self):
        geoms = self.selected_geometries()
        if not geoms:
            self.help_text.setPlainText("No geometry selected. Check one or more geometries above.")
            return
        # show help for first selected + a hint if multiple
        main = geoms[0]
        text = GEOM_HELP.get(main, "")
        if len(geoms) > 1:
            text += "\n\nAlso selected:\n- " + "\n- ".join(geoms[1:])
        text += (
            "\n\nGlobal parameters:\n"
            "- θ is in radians (if plotting vs theta)\n"
            "- q = 2k sin(θ/2), k = 2π n / λ\n"
            "- Units: use consistent length units for R, Rg, λ (e.g. nm)\n"
        )
        self.help_text.setPlainText(text)

    # ---- Parameter getters / x-grid

    def _get_params(self):
        R = float(self.radius.value())
        Rg = float(self.rg.value())
        dR = float(self.shell_thickness.value())
        lam = float(self.wavelength.value())
        nm = float(self.nm.value()) / 100.0
        return R, Rg, dR, lam, nm

    def _make_x(self, lam, nm):
        N = int(self.theta_n.value())

        if self.xmode.currentText().startswith("theta"):
            tmin = self.theta_min.value() / 100.0
            tmax = self.theta_max.value() / 100.0
            if tmax <= tmin:
                tmax = tmin + 0.01
            theta = np.linspace(tmin, tmax, N)
            q = q_from_theta(theta, nm, lam)
            return theta, q, "theta [rad]"
        else:
            qmin = self.q_min.value() / 1000.0
            qmax = self.q_max.value() / 1000.0
            if qmax <= qmin:
                qmax = qmin + 1e-3
            q = np.linspace(qmin, qmax, N)
            return q, q, "q [1/length]"

    # ---- Core plotting
    def _compute_geom(self, geom, q, R, Rg, dR):
        """
        Returns (F, I) where:
        - F: amplitude F(q) if available, else None
        - I: intensity |F(q)|^2 (or equivalent normalized form factor) if available, else None

        Note: Gaussian coil (Debye) and Guinier provide intensity I(q) only.
        """
        if geom == "Sphere":
            F = F_sphere_q(q, R)
            I = np.abs(F) ** 2
            return F, I

        if geom == "Thin spherical shell":
            F = F_thin_spherical_shell_q(q, R)
            I = np.abs(F) ** 2
            return F, I

        if geom == "Thick spherical shell":
            r_in = max(1e-9, R - dR)
            r_out = R
            F = F_thick_spherical_shell_q(q, r_in, r_out)
            I = np.abs(F) ** 2
            return F, I

        if geom == "Disk":
            F = F_disk_q(q, R)  # quick-look: q_perp ≈ q
            I = np.abs(F) ** 2
            return F, I

        if geom == "Gaussian coil (Debye)":
            return None, P_gaussian_coil_q(q, Rg)

        if geom == "Guinier":
            return None, P_guinier_q(q, Rg)

        return None, None

    def update_plot(self):
        geoms = self.selected_geometries()
        ax = self.canvas.ax
        ax.clear()

        if not geoms:
            ax.set_title("No geometry selected")
            self.canvas.draw()
            return

        R, Rg, dR, lam, nm = self._get_params()
        x, q, xlabel = self._make_x(lam, nm)

        ymode = self.ytype.currentText()

        plotted_any = False
        for geom in geoms:
            F, I = self._compute_geom(geom, q, R, Rg, dR)

            if ymode == "Re[F(q)]":
                if F is None:
                    continue
                y = np.real(F)

            elif ymode == "|F(q)|":
                if F is None:
                    continue
                y = np.abs(F)

            elif ymode == "|F(q)|²":
                if I is None:
                    continue
                y = np.asarray(I, dtype=float)

            else:
                continue

            ax.plot(x, y, label=geom, linewidth=1.75)
            plotted_any = True

        if not plotted_any:
            ax.set_title("No plottable geometry for selected Y mode")
            self.canvas.draw()
            return

        ax.set_xlabel(xlabel)
        ax.set_ylabel(ymode)  # ymode already equals the desired label text

        # x scale
        if self.xmode.currentText() == "q" and self.logx.isChecked():
            ax.set_xscale("log")
        else:
            ax.set_xscale("linear")

        # y scale (avoid log of negative values)
        if self.logy.isChecked() and ymode == "Re[F(q)]":
            ax.set_yscale("linear")
        else:
            ax.set_yscale("log" if self.logy.isChecked() else "linear")

        ax.set_title(f"λ={lam:g}, n={nm:.3g}, R={R:g}, Rg={Rg:g}, ΔR={dR:g}")
        ax.grid(True, alpha=0.25)
        ax.legend(loc="best", fontsize=10)
        self.canvas.draw()


    def save_plot(self):
        """
        Save the current matplotlib figure as a PNG.
        """
        filename, _ = QtWidgets.QFileDialog.getSaveFileName(
            self,
            "Save plot",
            "",
            "PNG files (*.png)"
        )

        if not filename:
            return

        if not filename.lower().endswith(".png"):
            filename += ".png"

        # High-quality export
        self.canvas.figure.savefig(
            filename,
            dpi=300,
            bbox_inches="tight"
        )

def run():
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    app.exec_()