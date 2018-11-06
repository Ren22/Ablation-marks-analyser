import sys, os, json, copy, re
from PyQt5.QtWidgets import *
from PyQt5.QtCore import Qt, pyqtSignal
from PyQt5.QtGui import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas, \
    NavigationToolbar2QT
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib.image as mpimg
import numpy as np
import matplotlib.colors
from matplotlib.gridspec import GridSpec

class State:
    def __init__(self):
        self.dataUploaded = False
        self.history = []
        self.log_scale = True

class LoadFiles(QWidget):

    updateDataSignal = pyqtSignal(dict)

    def __init__(self, state, canvas, mols_list):
        super().__init__()
        self.componentState = state
        self.canvas = canvas
        self.mols_list = mols_list

        #The code below can be potentially substituted by a button to open the file from another directory
        filename = '/home/renat/EMBL/spaceM_Luca/linux/testSamples/c2_SELECTED/Analysis/ili/sm_annotation_detections.csv'
        df = pd.read_csv(filename, delimiter=',')
        self.mols_df = df.iloc[:, 5:]
        self.mol_names = list(self.mols_df.columns.values)
        self.len_XY = int(np.sqrt(np.shape(df)[0]))
        self.arrX, self.arrY, self.randMol = df.X, df.Y, np.reshape(df[self.mol_names[0]], (self.len_XY, self.len_XY)).ravel()

        imagepath = '/home/renat/EMBL/spaceM_Luca/linux/testSamples/c2_SELECTED/Analysis/ili/FLUO_crop_bin1x1.png'
        self.imgplt = mpimg.imread(imagepath)

        cellProfImgPath = '/home/renat/EMBL/spaceM_Luca/linux/testSamples/c2_SELECTED/Analysis/scAnalysis/Molecular_features/marks_flitered_fluo.npy'
        self.cellProfImg = np.load(cellProfImgPath)
        self.get13vals = self.cellProfImg[13]

        self.componentState.dataUploaded = True
        self.mols_list.update_mols_df(self.mol_names[0], self.mols_df)
        self.canvas.arrX, self.canvas.arrY, self.canvas.arrZ = self.arrX, self.arrY, self.randMol
        self.canvas.img = self.imgplt
        self.canvas.val13 = self.get13vals
        self.canvas.clean_n_plot(arrX=self.arrX, arrY=self.arrY, arrZ=self.randMol, img=self.imgplt, val13=self.get13vals)

    def __call__(self):
        return {
            'arrX': self.arrX,
            'arrY': self.arrY,
            'randMol': self.randMol,
            'imgplt': self.imgplt,
            'mols_df': self.mols_df,
        }

    def view(self):
        # Add a file button
        self.openBtn = QPushButton('Open file', self)
        self.openBtn.move(10, 10)
        self.openBtn.adjustSize()
        self.openBtn.clicked.connect(lambda: self.openFile())

class NavigationToolbar(NavigationToolbar2QT):
    toolitems = [t for t in NavigationToolbar2QT.toolitems if
                 t[0] in ('Home', 'Back', 'Forward', 'Pan', 'Zoom', 'Subplots')]

class Mols_Section(QWidget):

    signal = pyqtSignal(str, list)

    def __init__(self, state):
        super().__init__()
        self.componentState = state
        self.cb = QComboBox(self)
        self.mols_df = []
        self.molecule = self.cb.currentText()
        self.cb.activated[str].connect(self.handle_molName)

    def update_mols_df(self, molecule, mols_df_new):
        self.molecule = molecule
        self.mols_df = mols_df_new
        if self.componentState:
            self.cb.addItems(self.mols_df.columns.values)

    def handle_molName(self, mol_name):
        self.signal.emit(mol_name, self.mols_df[mol_name].tolist())

    @property
    def currMolecule(self):
        return self.molecule

    @property
    def currMoleculeValues(self):
        if self.componentState:
            return self.mols_df[self.molecule]
        else:
            return []

class Tabs(QWidget):
    flippingSignal = pyqtSignal(str)
    spotsSizeSignal = pyqtSignal(int)
    degreeSignal = pyqtSignal(int)
    logScaleSignal = pyqtSignal(bool)
    exportConfigsSingal = pyqtSignal()

    def __init__(self, state, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)
        self.componentState = state

        self.tabs = QTabWidget()
        self.tab1 = QWidget()

        # GroupBox
        self.groupBox1 = QGroupBox("Rotate clockwise: (Rotation is always relative to the initial position)")
        self.radio90 = QRadioButton("90 deg")
        self.radio180 = QRadioButton("180 deg")
        self.radio270 = QRadioButton("270 deg")
        self.vbox = QVBoxLayout()
        self.vbox.addWidget(self.radio90)
        self.vbox.addWidget(self.radio180)
        self.vbox.addWidget(self.radio270)
        self.groupBox1.setLayout(self.vbox)
        self.radio90.toggled.connect(lambda: self.gbDegSignal(90))
        self.radio180.toggled.connect(lambda: self.gbDegSignal(180))
        self.radio270.toggled.connect(lambda: self.gbDegSignal(270))

        # Checkboxes
        self.flippingLabel = QLabel(self)
        self.flippingLabel.setText('Flip the image:')
        self.cb_flipud = QCheckBox("Flip ion image OX", self)
        self.cb_fliplr = QCheckBox("Flip ion image OY", self)
        self.cb_transp = QCheckBox("Transponse ion image", self)
        self.logscale = QCheckBox("Logarithmic scale")
        self.logscale.setChecked(True)

        # Checkboxes & Buttons handlers
        self.cb_flipud.clicked.connect(self.flipUDsignal)
        self.cb_fliplr.clicked.connect(self.flipLRsignal)
        self.cb_transp.clicked.connect(self.transpSignal)
        self.logscale.stateChanged.connect(self.stateC)

        # Sliders
        self.sliderLabel = QLabel(self)
        self.sliderLabel.setText('Choose spots size:')
        self.slider = QSlider(Qt.Horizontal)
        self.slider.setFocusPolicy(Qt.StrongFocus)
        self.slider.setTickPosition(QSlider.TicksBothSides)
        self.slider.setTickInterval(10)
        self.slider.setSingleStep(10)
        self.slider.setValue(60)
        self.slider.setMinimum(0)
        self.slider.setMaximum(120)

        # Buttons
        self.btn_reset = QPushButton('Reset ion image')
        self.btn_reset.clicked.connect(self.reset_ion_image)
        self.btn_save_configs = QPushButton('Save configurations')
        self.btn_save_configs.clicked.connect(self.export_configs)

        #Slider handlers
        self.slider.valueChanged.connect(self.spotsSizeSignal)

        # Add tabs
        self.tabs.addTab(self.tab1, "Transformations")
        # self.tabs.addTab(self.tab3, "Mapping")

        # Create 1st tab
        self.tab1.layout = QVBoxLayout(self)
        self.tab1.layout.addWidget(self.groupBox1)
        self.tab1.layout.addWidget(self.flippingLabel)
        self.tab1.layout.addWidget(self.cb_flipud)
        self.tab1.layout.addWidget(self.cb_fliplr)
        self.tab1.layout.addWidget(self.cb_transp)
        self.tab1.layout.addWidget(self.btn_reset)
        self.tab1.layout.addWidget(self.sliderLabel)
        self.tab1.layout.addWidget(self.slider)
        self.tab1.layout.addWidget(self.logscale)
        self.tab1.layout.addWidget(self.btn_save_configs)

        self.tab1.layout.addStretch()
        self.tab1.setLayout(self.tab1.layout)

        self.layout.addWidget(self.tabs)

    def flipLRsignal(self):
        if self.componentState.dataUploaded:
            if 'fliplr' in self.componentState.history:
                self.componentState.history.remove('fliplr')
                self.componentState.history = self.remove_double_int(self.componentState.history)
            else:
                self.componentState.history.append('fliplr')
            self.flippingSignal.emit('fliplr')

            print(self.componentState.history)

    def flipUDsignal(self):
        if self.componentState.dataUploaded:
            if 'flipud' in self.componentState.history:
                self.componentState.history.remove('flipud')
                self.componentState.history = self.remove_double_int(self.componentState.history)
            else:
                self.componentState.history.append('flipud')
            self.flippingSignal.emit('flipud')

            print(self.componentState.history)

    def transpSignal(self):
        if self.componentState.dataUploaded:
            if 'transp' in self.componentState.history:
                self.componentState.history.remove('transp')
                self.componentState.history = self.remove_double_int(self.componentState.history)
            else:
                self.componentState.history.append('transp')
            self.flippingSignal.emit('transp')

            print(self.componentState.history)

    def reset_ion_image(self):
        if self.componentState.dataUploaded:
            k = None
            self.cb_flipud.setChecked(False)
            self.cb_fliplr.setChecked(False)
            self.cb_transp.setChecked(False)
            self.radio90.setAutoExclusive(False)
            self.radio180.setAutoExclusive(False)
            self.radio270.setAutoExclusive(False)
            self.radio90.setChecked(False)
            self.radio180.setChecked(False)
            self.radio270.setChecked(False)
            self.radio90.setAutoExclusive(True)
            self.radio180.setAutoExclusive(True)
            self.radio270.setAutoExclusive(True)
            self.flippingSignal.emit('reset')
            self.componentState.history.clear()
            print(self.componentState.history)

    def spotsSize(self, sval):
        if self.componentState.dataUploaded:
            self.spotsSizeSignal.emit(sval)

    def gbDegSignal(self, deg):
        # TODO: there's a bug below: each time rotation angle is changed the previous angle is double-
        # called because of toggled radio-button value
        deg_to_k = {
            90: 1,
            180: 2,
            270: 3
        }
        if self.componentState.dataUploaded:
            if self.componentState.history and type(self.componentState.history[-1]) is int:
                self.componentState.history.pop()
                self.componentState.history.append(deg_to_k[deg])
            else:
                self.componentState.history.append(deg_to_k[deg])
            print(self.componentState.history)
            self.flippingSignal.emit('flip_clwise')

    def stateC(self):
        if self.componentState.dataUploaded:
            if self.logscale.isChecked():
                self.componentState.log_scale = True
                self.logScaleSignal.emit(True)
            else:
                self.componentState.log_scale = False
                self.logScaleSignal.emit(False)

    def export_configs(self):
        print(self.componentState.history)
        with open("AM_analysis_config.txt", "w") as f:
            f.write(json.dumps(self.componentState.history))

    @staticmethod
    def remove_double_int(history):
        i = 0
        if len(history) == 2 and type(history[i]) is int and type(history[i + 1]) is int:
            history.remove(history[i])
        while len(history) > 2 and i < len(history) - 2:
            if type(history[i]) is int and type(history[i + 1]) is int:
                history.remove(history[i])
            i = i + 1
        return history

class Window(QMainWindow):

    def __init__(self):
        super().__init__()
        self.state = State()
        self.central_widget = QWidget(self) #central widget needs to be created first
        self.setCentralWidget(self.central_widget)
        self.createOtherWidgets(parent=self.central_widget)
        self.layout = View(self.getWidgets)
        self.central_widget.setLayout(self.layout.lt)
        self.methods = Controller(self.getWidgets, self.state)
        self.initWinUI()
        np.set_printoptions(threshold=np.nan)

    def initWinUI(self):
        self.setWindowTitle('Ablation marks registration analysis')
        self.setGeometry(100, 100, 1200, 800)
        self.show()

    def createOtherWidgets(self, parent):
        self.canvas = MatplotlibArea(self.state)
        self.navigation = NavigationToolbar(self.canvas, parent)
        self.tabs = Tabs(self.state, parent)
        self.mols = Mols_Section(self.state)
        self.data = LoadFiles(self.state, self.canvas, self.mols)

    @property
    def getWidgets(self):
        return {
            'canvas': self.canvas,
            'navigation': self.navigation,
            'tabs': self.tabs,
            'mols': self.mols,
            'data': self.data
        }

class Controller:

    def __init__(self, widgets=None, state=None):
        self.componentState = state
        self.canvas = widgets['canvas']
        self.navigation = widgets['navigation']
        self.tabs = widgets['tabs']
        self.data = widgets['data']
        self.mols = widgets['mols']
        self.spotVal = self.canvas.spotsize
        self.normPlot = self.canvas.norm
        self.currMol = self.mols.currMolecule
        self.currMolVals = self.mols.currMoleculeValues
        self.init_currMolVals = copy.deepcopy(self.currMolVals)
        self.rf = int(self.currMolVals.shape[0]**0.5) #rf = reshaping factor

        # Signal handlers
        self.mols.signal.connect(self.updateCanvas)
        self.tabs.flippingSignal.connect(self.flip)
        self.tabs.spotsSizeSignal.connect(self.changeSpotsSize)
        self.tabs.logScaleSignal.connect(self.enableLogScale)
        self.k = 0

        # Define transforming functions
        self.trans_func = {}
        self.trans_func['fliplr'] = np.fliplr
        self.trans_func['flipud'] = np.flipud
        self.trans_func['transp'] = np.transpose
        self.trans_func['flip_clwise'] = np.rot90

    def apply_transforms(self):
        curr_mol_vals = self.init_currMolVals
        for i in self.componentState.history:
            if type(i) is int:
                curr_mol_vals = self.trans_func['flip_clwise'](np.array(curr_mol_vals).reshape(self.rf, self.rf), k=i).ravel()
            else:
                curr_mol_vals = self.trans_func[i](np.array(curr_mol_vals).reshape(self.rf, self.rf)).ravel()
        return curr_mol_vals

    def redraw(self):
        self.canvas.clean_n_plot(arrX=self.canvas.arrX,
                                 arrY=self.canvas.arrY,
                                 arrZ=self.currMolVals,
                                 s=self.spotVal,
                                 val13=self.canvas.val13,
                                 img=self.canvas.img)

    def updateCanvas(self, mol_name, mols_vals):
        self.currMol = mol_name
        self.currMolVals = np.reshape(mols_vals, (self.rf, self.rf)).ravel()
        self.init_currMolVals = copy.deepcopy(self.currMolVals)
        if self.componentState.history:
            self.currMolVals = self.apply_transforms()
        self.canvas.norm = matplotlib.colors.LogNorm() if self.componentState.log_scale is True else None
        self.redraw()

    def flip(self, flipside):

        if flipside == 'fliplr':
            self.currMolVals = self.apply_transforms()
            self.redraw()

        if flipside == 'flipud':
            self.currMolVals = self.apply_transforms()
            self.redraw()

        if flipside == 'flip_clwise':
            self.currMolVals = self.apply_transforms()
            self.redraw()

        if flipside == 'transp':
            self.currMolVals = self.apply_transforms()
            self.redraw()

        if flipside == 'reset':
            self.currMolVals = self.init_currMolVals
            self.redraw()
            self.componentState.flip_horizontally = False
            self.componentState.flip_vertically = False
            self.componentState.transposed = False

    def changeSpotsSize(self, sval):
        if self.componentState.dataUploaded:
            self.spotVal = sval
            self.redraw()

    def enableLogScale(self, val):
        self.canvas.norm = matplotlib.colors.LogNorm() if self.componentState.log_scale is True else None
        if val and self.componentState.dataUploaded:
            self.redraw()
        else:
            self.redraw()

class View:

    def __init__(self, widgets):
        self.layout = QGridLayout()
        self.canvas = widgets['canvas']
        self.navigation = widgets['navigation']
        self.tabs = widgets['tabs']
        self.mols = widgets['mols']
        self.data = widgets['data']

        self.layout.setColumnStretch(0, 6)
        self.layout.setColumnStretch(1, 6)
        self.layout.addWidget(self.navigation, 0, 0)
        self.layout.addWidget(self.canvas, 1, 0)
        self.layout.addWidget(self.tabs, 1, 1)
        self.layout.addWidget(self.mols, 0, 1)

    @property
    def lt(self):
        return self.layout

class MatplotlibArea(FigureCanvas):

    def __init__(self, state, width=5, height=4, dpi=100):
        self.componentState = state
        self.arrX, self.arrY, self.arrZ = [], [], []
        self.img = None
        self.val13 = None
        self.spotsize = 20
        self.norm = matplotlib.colors.LogNorm() if self.componentState.log_scale is True else None
        self.fig = plt.figure(figsize=(width, height), dpi=dpi)
        self.limX = ()
        self.limY = ()

        gs = GridSpec(3, 2)
        self.ax1 = plt.subplot(gs[:-1, :])
        self.ax2 = plt.subplot(gs[-1, 0])
        self.ax3 = plt.subplot(gs[-1, 1])
        super(MatplotlibArea, self).__init__(self.fig)

    def plot(self, arrX, arrY, arrZ, img, val13, s):
        if img is not None:
            self.ax1.imshow(img)
        self.ax1.scatter(arrX, arrY, s, c=arrZ, norm=self.norm, edgecolor='')
        self.ax1.callbacks.connect('xlim_changed', self.on_xlims_change)
        self.ax1.callbacks.connect('ylim_changed', self.on_ylims_change)
        rf = int(val13.shape[0]**0.5)
        if val13 is not None:
            self.ax2.imshow(np.reshape(val13, (rf, rf)))
        # self.ax2.axis('off')
        # self.ax3.axis('off')
        self.ax3.imshow(np.reshape(arrZ, (rf, rf)), norm=self.norm)

    def on_xlims_change(self, axes):
        self.limX = axes.get_xlim()

    def on_ylims_change(self, axes):
        self.limY = axes.get_ylim()

    def clean_n_plot(self, arrX, arrY, arrZ, img, val13, s=None):
        self.ax1.cla()
        self.ax2.cla()
        self.ax3.cla()
        if self.limX and self.limY:
            self.ax1.set_xlim(self.limX)
            self.ax1.set_ylim(self.limY)
        self.plot(arrX, arrY, arrZ, img, val13, s)
        # use canvas draw rather than artist draw, issue https://sourceforge.net/p/matplotlib/mailman/message/23209300/
        self.ax1.figure.canvas.draw()
        self.ax2.figure.canvas.draw()
        self.ax3.figure.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = Window()
    sys.exit(app.exec_())