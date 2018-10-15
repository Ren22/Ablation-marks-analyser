import sys, os, json
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
        self.flip_horizontally = False
        self.flip_vertically = False
        self.rotation_deg = 0
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
    flippingSignal = pyqtSignal(str, int)
    spotsSizeSignal = pyqtSignal(int)
    degreeSignal = pyqtSignal(int)
    logScaleSignal = pyqtSignal(bool)

    def __init__(self, state, parent):
        super(QWidget, self).__init__(parent)
        self.layout = QVBoxLayout(self)
        self.componentState = state
        self.deg_old = -1

        self.tabs = QTabWidget()
        self.tab1 = QWidget()
        self.tab2 = QWidget()
        self.tab3 = QWidget()

        # GroupBox
        self.groupBox1 = QGroupBox("Rotate clockwise:")
        self.radio0 = QRadioButton("0 deg")
        self.radio90 = QRadioButton("90 deg")
        self.radio180 = QRadioButton("180 deg")
        self.radio270 = QRadioButton("270 deg")
        self.vbox = QVBoxLayout()
        self.vbox.addWidget(self.radio0)
        self.vbox.addWidget(self.radio90)
        self.vbox.addWidget(self.radio180)
        self.vbox.addWidget(self.radio270)
        self.groupBox1.setLayout(self.vbox)
        self.radio0.toggled.connect(lambda: self.gbDegSignal(0))
        self.radio90.toggled.connect(lambda: self.gbDegSignal(90))
        self.radio180.toggled.connect(lambda: self.gbDegSignal(180))
        self.radio270.toggled.connect(lambda: self.gbDegSignal(270))

        # Checkboxes
        self.flippingLabel = QLabel(self)
        self.flippingLabel.setText('Flip the image: \nWarning: this won\'t be preserved when rotating image')
        self.cb_flipud = QCheckBox("Flip image OX", self)
        self.cb_fliplr = QCheckBox("Flip image OY", self)
        self.logscale = QCheckBox("Logarithmic scale")
        self.logscale.setChecked(True)
        # self.state.log_scale

        # Checkboxes & Buttons handlers
        self.cb_flipud.clicked.connect(self.flipUDsignal)
        self.cb_fliplr.clicked.connect(self.flipLRsignal)
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
        self.btn_save_configs = QPushButton('Save configurations')
        self.btn_save_configs.clicked.connect(self.export_configs)

        #Slider handlers
        self.slider.valueChanged.connect(self.spotsSizeSignal)

        # Add tabs
        self.tabs.addTab(self.tab1, "Transformations")
        self.tabs.addTab(self.tab3, "Mapping")

        # Create 1st tab
        self.tab1.layout = QVBoxLayout(self)
        self.tab1.layout.addWidget(self.flippingLabel)
        self.tab1.layout.addWidget(self.cb_flipud)
        self.tab1.layout.addWidget(self.cb_fliplr)
        self.tab1.layout.addWidget(self.groupBox1)
        # Copypasted for presentation only
        self.tab1.layout.addWidget(self.sliderLabel)
        self.tab1.layout.addWidget(self.slider)
        # Copypasted for presentation only
        self.tab1.layout.addWidget(self.logscale)
        self.tab1.layout.addWidget(self.btn_save_configs)


        self.tab1.layout.addStretch()
        self.tab1.setLayout(self.tab1.layout)

        # Create 2nd tab
        # self.tab2.layout = QVBoxLayout(self)
        # self.tab2.layout.addWidget(self.slider)
        # self.tab2.layout.addStretch()
        # self.tab2.setLayout(self.tab2.layout)

        self.layout.addWidget(self.tabs)
        # self.setLayout(self.layout)

    def flipLRsignal(self):
        if self.componentState.dataUploaded:
            k = 1
            self.flippingSignal.emit('fliplr', k)

    def flipUDsignal(self):
        if self.componentState.dataUploaded:
            k = 1
            self.flippingSignal.emit('flipud', k)

    def spotsSize(self, sval):
        if self.componentState.dataUploaded:
            self.spotsSizeSignal.emit(sval)

    def gbDegSignal(self, deg):
        # TODO: there's a bug below: each time rotation angle is changed the previous angle is double-
        # called because of toggled radio-button value
        if self.componentState.dataUploaded:
            if deg == 0 and deg != self.deg_old:
                k = 0
            elif deg == 90 and deg != self.deg_old:
                k = 1
            elif deg == 180 and deg != self.deg_old:
                k = 2
            elif deg == 270 and deg != self.deg_old:
                k = 3
            if deg != self.deg_old:
                print(deg, k)
                self.flippingSignal.emit('flip_clockwise', k)
            self.deg_old = deg

    def stateC(self):
        if self.componentState.dataUploaded:
            if self.logscale.isChecked():
                self.logScaleSignal.emit(True)
            else:
                self.logScaleSignal.emit(False)

    def export_configs(self):
        # if self.componentState['rotation_deg'] > 0:
        #     k = (self.componentState['rotation_deg']/90) % 4
        # if self.componentState['flip_horizontally'] is True:
        #     return json.dumps(self.componentState)
        with open("AM_analysis_config.txt", "w") as f:
            f.write(json.dumps(vars(self.componentState)))

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
        # if self.state.log_scale is True:
        self.normPlot = self.canvas.norm
        self.currMol = self.mols.currMolecule
        self.currMolVals = self.mols.currMoleculeValues
        self.rf = int(self.currMolVals.shape[0]**0.5) #rf = reshaping factor

        self.mols.signal.connect(self.updateCanvas)
        self.tabs.flippingSignal.connect(self.flip)
        self.tabs.spotsSizeSignal.connect(self.changeSpotsSize)
        self.tabs.logScaleSignal.connect(self.enableLogScale)
        self.k_to_full_rot = 0

    def updateCanvas(self, mol_name, mols_vals):
        self.currMol = mol_name
        self.lenXY = int(np.sqrt(np.shape(mols_vals)[0]))
        self.currMolVals = np.reshape(mols_vals, (self.lenXY, self.lenXY)).ravel()
        self.canvas.clean_n_plot(arrX=self.canvas.arrX,
                                 arrY=self.canvas.arrY,
                                 arrZ=self.currMolVals,
                                 s=self.spotVal,
                                 norm=self.normPlot,
                                 val13=self.canvas.val13,
                                 img=self.canvas.img)


    def flip(self, flipside, k=None):

        if flipside == 'fliplr' and self.componentState.dataUploaded:
            self.currMolVals = np.fliplr(np.array(self.currMolVals).reshape(self.rf, self.rf)).ravel()
            self.canvas.clean_n_plot(arrX=self.canvas.arrX,
                                     arrY=self.canvas.arrY,
                                     arrZ=self.currMolVals,
                                     s=self.spotVal,
                                     norm=self.normPlot,
                                     val13=self.canvas.val13,
                                     img=self.canvas.img)
            self.componentState.flip_vertically = not self.componentState.flip_vertically
            print(vars(self.componentState))

        if flipside == 'flipud' and self.componentState.dataUploaded:
            self.currMolVals = np.flipud(np.array(self.currMolVals).reshape(self.rf, self.rf)).ravel()
            self.canvas.clean_n_plot(arrX=self.canvas.arrX,
                                     arrY=self.canvas.arrY,
                                     arrZ=self.currMolVals,
                                     s=self.spotVal,
                                     norm=self.normPlot,
                                     val13=self.canvas.val13,
                                     img=self.canvas.img)
            self.componentState.flip_horizontally = not self.componentState.flip_horizontally
            print(vars(self.componentState))

        if flipside == 'flip_clockwise' and self.componentState.dataUploaded:
            self.currMolVals = np.rot90(np.array(self.currMolVals).reshape(self.rf, self.rf), k).ravel()
            self.canvas.clean_n_plot(arrX=self.canvas.arrX,
                                     arrY=self.canvas.arrY,
                                     arrZ=self.currMolVals,
                                     s=self.spotVal,
                                     norm=self.normPlot,
                                     val13=self.canvas.val13,
                                     img=self.canvas.img)
            # The trick below is to rotate image to the initial state (4-k) so the next run it will start from 0deg
            self.currMolVals = np.rot90(np.array(self.currMolVals).reshape(self.rf, self.rf), 4 - k).ravel()
            self.componentState.rotation_deg = 90*k
            print(vars(self.componentState))

    def changeSpotsSize(self, sval):
        if self.componentState.dataUploaded:
            self.spotVal = sval
            self.canvas.clean_n_plot(arrX=self.canvas.arrX,
                                     arrY=self.canvas.arrY,
                                     arrZ=self.currMolVals,
                                     s=self.spotVal,
                                     norm=self.normPlot,
                                     val13=self.canvas.val13,
                                     img=self.canvas.img)

    def enableLogScale(self, val):
        if val and self.componentState.dataUploaded:
            self.normPlot = matplotlib.colors.LogNorm()
            self.canvas.clean_n_plot(arrX=self.canvas.arrX,
                                     arrY=self.canvas.arrY,
                                     arrZ=self.currMolVals,
                                     s=self.spotVal,
                                     norm=self.normPlot,
                                     val13=self.canvas.val13,
                                     img=self.canvas.img,
                                     )
        else:
            self.normPlot = None
            self.canvas.clean_n_plot(arrX=self.canvas.arrX,
                                     arrY=self.canvas.arrY,
                                     arrZ=self.currMolVals,
                                     s=self.spotVal,
                                     norm=self.normPlot,
                                     val13=self.canvas.val13,
                                     img=self.canvas.img,
                                     )

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
        # self.layout.addWidget(self.data, 0, 2)
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
        self.norm = None
        self.fig = plt.figure(figsize=(width, height), dpi=dpi)
        self.limX = ()
        self.limY = ()

        gs = GridSpec(3, 2)
        self.ax1 = plt.subplot(gs[:-1, :])
        self.ax2 = plt.subplot(gs[-1, 0])
        self.ax3 = plt.subplot(gs[-1, 1])
        super(MatplotlibArea, self).__init__(self.fig)

    def plot(self, arrX, arrY, arrZ, img, val13, s, norm):
        if img is not None:
            self.ax1.imshow(img)
        self.ax1.scatter(arrX, arrY, s, c=arrZ, norm=norm, edgecolor='')
        self.ax1.callbacks.connect('xlim_changed', self.on_xlims_change)
        self.ax1.callbacks.connect('ylim_changed', self.on_ylims_change)
        rf = int(val13.shape[0]**0.5)
        if val13 is not None:
            self.ax2.imshow(np.reshape(val13, (rf, rf)).T)
            # self.ax2.axis('off')
        # self.ax3.axis('off')
        self.ax3.imshow(np.reshape(arrZ, (rf, rf)))

    def on_xlims_change(self, axes):
        self.limX = axes.get_xlim()

    def on_ylims_change(self, axes):
        self.limY = axes.get_ylim()

    def clean_n_plot(self, arrX, arrY, arrZ, img, val13, norm=None, s=None):
        self.ax1.cla()
        self.ax2.cla()
        self.ax3.cla()
        if self.limX and self.limY:
            self.ax1.set_xlim(self.limX)
            self.ax1.set_ylim(self.limY)
        self.plot(arrX, arrY, arrZ, img, val13, s, norm)
        # if s is None:
        #     self.plot(arrX, arrY, arrZ, img)
        # else:
        #
        # use canvas draw rather than artist draw, issue https://sourceforge.net/p/matplotlib/mailman/message/23209300/
        self.ax1.figure.canvas.draw()
        self.ax2.figure.canvas.draw()
        self.ax3.figure.canvas.draw()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    main = Window()
    sys.exit(app.exec_())