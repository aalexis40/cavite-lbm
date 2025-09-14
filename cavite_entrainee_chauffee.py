import sys
import numpy as np
from PyQt5.QtWidgets import QApplication, QMainWindow, QLabel, QVBoxLayout, QPushButton, QWidget, QLineEdit, QSizePolicy, QSlider, QProgressBar, QFileDialog
from PyQt5.QtGui import QIcon, QPixmap
from PyQt5.QtCore import QTimer, Qt, pyqtSignal, QObject
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.animation import FuncAnimation
from abc import ABC, abstractmethod
import copy
import matplotlib.animation as animation

class InterfaceGraphique(QMainWindow): # fenêtre graphique principale
    def __init__(self, m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D):
        super().__init__()
        self.setWindowTitle("Cavité enteainée chauffée")
        self.sub_interface = None  
        self.m = m  # taille domaine
        self.delta = delta  #pas de grille
        self.tau = tau  # paramètre de relaxation
        self.rho0 = rho0  # densité initiale
        self.chi0 = chi0  # chi initial 
        self.v0 = v0  # vitesse initiale
        self.U0 = U0  # vitesse du bord supérieur
        self.nu = nu  # Viscosité 
        self.Re = Re  # nombre de Reynolds
        self.Pr = Pr  # nombre de Prandtl
        self.D = D  # coeff de diffusion
        
        central_widget = QWidget()
        self.setCentralWidget(central_widget)
        layout = QVBoxLayout()
        central_widget.setLayout(layout)

        self.creer_bouton(layout, "Visualiser animation")
        self.creer_bouton(layout, "Visualiser courbes isothermes")
        self.creer_bouton(layout, "Visualiser streamlines")

        self.bouton_quitter = QPushButton("Quitter", self)
        self.bouton_quitter.clicked.connect(self.close)
        layout.addWidget(self.bouton_quitter)

        with open('style.css', 'r') as f:
            self.setStyleSheet(f.read())
            
    def creer_bouton(self, layout, button_text):# méthode pour créer un bouton avec une entrée
        frame = QWidget()
        frame_layout = QVBoxLayout()
        frame.setLayout(frame_layout)

        bouton = QPushButton(button_text)
        bouton.clicked.connect(lambda: self.open_sub_interface(button_text)) # ouvrir la fenêtre graphique secondaire lors de l'appui sur le bouton
        frame_layout.addWidget(bouton)

        layout.addWidget(frame)

    def open_sub_interface(self, button_text): # méthode pour ouvrir la fenêtre secondaire
        self.sub_interface = SubInterface(self, button_text)
        self.sub_interface.show()

class SubInterface(QWidget): # fenêtre graphique secondaire pour l'animation
    def __init__(self, parent, button_text):
        super().__init__()
        self.setWindowTitle(f"Options pour {button_text}")
        self.parent = parent
        self.button_text = button_text


        self.layoutSub = QVBoxLayout()

        self.fig, self.ax = plt.subplots()
        self.ax.set_axis_off()
        self.canvas = FigureCanvas(self.fig)
        self.layoutSub.addWidget(self.canvas)

        
        self.label = QLabel("Vitesse x 1")
        self.layoutSub.addWidget(self.label)


        self.slider = QSlider(Qt.Horizontal)
        self.slider.setMinimum(1)
        self.slider.setMaximum(20)
        self.slider.setValue(1)
        self.slider.setTickInterval(20)
        self.slider.setTickPosition(QSlider.TicksBelow)
        self.slider.valueChanged.connect(self.update_value)
        self.layoutSub.addWidget(self.slider)

        self.setLayout(self.layoutSub)

        label = QLabel(f"Options pour {button_text}")
        self.layoutSub.addWidget(label)

        nb_steps_label = QLabel("Nombre d'itérations:")
        self.layoutSub.addWidget(nb_steps_label)

        self.nb_steps_entry = QLineEdit()
        self.layoutSub.addWidget(self.nb_steps_entry)


        self.bouton_afficher = QPushButton("Afficher")
        self.bouton_afficher.clicked.connect(self.afficher)
        self.layoutSub.addWidget(self.bouton_afficher)

        self.bouton_enregistrer = QPushButton("Enregistrer")
        self.bouton_enregistrer.clicked.connect(self.enregistrer)
        self.layoutSub.addWidget(self.bouton_enregistrer)
        self.bouton_enregistrer.hide()
        
        self.bouton_quitter = QPushButton("Fermer", self)
        self.bouton_quitter.clicked.connect(self.close)
        self.layoutSub.addWidget(self.bouton_quitter)

        self.chi_results = None
        self.bloc = False
        self.crb = False
        self.animation_terminee = False
        self.chi_results, self.vitesses_results, self.rho_results = None, None, None
        self.steps = None

        with open('style.css', 'r') as f:
            self.setStyleSheet(f.read())

        self.progress_bar = QProgressBar()
        self.layoutSub.addWidget(self.progress_bar)


        self.progress_bar.setRange(0, 100)

    def update_value(self, value):
        self.label.setText(f"Vitesse x {value}")

    def afficher(self):
        if not self.bloc:
            self.animation_terminee = False
            self.bloc = not self.bloc

            self.bouton_afficher.setText("Stop")

            self.steps = int(self.nb_steps_entry.text())

            self.progress_bar.setValue(0)  
            self.progress_bar.setMaximum(self.steps) 

            simu = Simulation(self.steps, self.parent.m, self.parent.delta, self.parent.tau, 
                              self.parent.rho0, self.parent.chi0, self.parent.v0, self.parent.U0, 
                              self.parent.nu, self.parent.Re, self.parent.Pr, self.parent.D)
            
            simu.set_nb_iteration(self.steps)
            simu.set_vitesse(self.slider.value())
            simu.set_intG(self)
            self.chi_results, self.vitesses_results, self.rho_results = simu.simul()

            # réinitialiser la figure et les axes
            self.fig.clear()
            self.ax = self.fig.add_subplot(111)
            self.ax.set_axis_off()

            # self.current_frame = 0
            # self.animation_timer = QTimer()
            if self.button_text == "Visualiser animation":
                self.crb = False
                self.animation = FuncAnimation(self.fig, self.update_courbes_iso, frames=len(self.chi_results), interval=0.2, repeat=False)

            elif self.button_text == "Visualiser courbes isothermes": # animation avec QTimer
                self.crb = True
                # self.animation_timer.timeout.connect(self.update_animation_frame) 
                # self.animation_timer.start(20) 
                self.animation = FuncAnimation(self.fig, self.update_courbes_iso, frames=len(self.chi_results), interval=0.2, repeat=False)
                
            elif self.button_text == "Visualiser streamlines": # animation avec FuncAnimation
                self.animation = FuncAnimation(self.fig, self.update_streamlines, frames=len(self.chi_results), interval=0.2, repeat=False)
                
            # Cacher les boutons
            self.label.hide()
            self.slider.hide()
            self.nb_steps_entry.hide()
        else:
            self.bloc = False
            if not self.animation_terminee:
                self.animation.event_source.stop()
            self.progress_bar.setValue(0)  
            self.bouton_afficher.setText("Afficher")
            self.show_widgets()



    # def update_animation_frame(self): # animation QTimer pour les courbes isothermes
    #     self.update_courbes_iso(self.current_frame)
    #     self.current_frame += 1
    #     if self.current_frame >= len(self.chi_results):
    #         self.animation_timer.stop()  

    def update_courbes_iso(self, frame): # animation pour les courbes isothermes
        self.ax.clear()
        self.ax.imshow(self.chi_results[frame], cmap='hot')
        if self.crb:
            self.ax.contour(self.chi_results[frame], colors='black', linewidths=0.5, levels=[0.2,0.25,0.3,0.4,0.5,0.6,0.7])
        self.ax.set_axis_off()
        self.canvas.draw()
        if frame == len(self.chi_results)-1:
            self.animation_terminee = True

    def update_streamlines(self, frame): # animation pour les streamlines

        self.ax.clear()
        self.ax.imshow(self.chi_results[frame], cmap='hot')
        self.ax.streamplot(np.arange(self.vitesses_results[frame][0].shape[1]), np.arange(self.vitesses_results[frame][1].shape[0]), self.vitesses_results[frame][0]*self.rho_results[frame], self.vitesses_results[frame][1]*(-self.rho_results[frame]), color='blue', density=0.5)
        self.ax.set_axis_off()
        self.canvas.draw()
        if frame == len(self.chi_results)-1:
            self.animation_terminee = True

    def show_widgets(self): # afficher les widgets cachés lors de l'animation
        # Afficher les boutons
        self.label.show()
        self.slider.show()
        self.nb_steps_entry.show()
        self.bouton_enregistrer.show()

    def enregistrer(self):
        options = QFileDialog.Options()
        file_name, _ = QFileDialog.getSaveFileName(self, "Enregistrer l'animation", "", "Video Files (*.mp4)", options=options)
        if file_name:
            DPI=90
            writer = animation.FFMpegWriter(fps=30, bitrate=5000)
            self.animation.save(file_name, writer=writer, dpi=DPI)
            self.bouton_enregistrer.hide()

class LBM:
    def __init__(self, m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D):
        # Initialisation des paramètres du modèle LBM
        self.m = m  # taille domaine
        self.delta = delta  #pas de grille
        self.tau = tau  # paramètre de relaxation
        self.rho0 = rho0  # densité initiale
        self.chi0 = chi0  # chi initial 
        self.v0 = v0  # vitesse initiale
        self.U0 = U0  # vitesse du bord supérieur
        self.nu = nu  # Viscosité 
        self.omega_m = 1 / (3 * nu + 0.5)  # paramètre de relaxation pour la vitesse
        self.Re = Re  # nombre de Reynolds
        self.Pr = Pr  # nombre de Prandtl
        self.D = D  # coeff de diffusion
        self.omega_s = 1 / (3 * D + 0.5)  # paramètre de relaxation pour la température

        self.chi = np.zeros((m, m))  # température
        self.rho = np.ones((m, m)) * rho0  # densité
        self.v = np.zeros((2,m, m))  # vitesse x

        # directions de la vitesse D2Q9
        self.c = np.array([[0, 0], [1, 0], [0, 1], [-1, 0], [0, -1], [1, 1], [-1, 1], [-1, -1], [1, -1]])
        self.c_s2 = 1/3  # (vitesse du son)²

        # Poids directions vitesse
        self.t = np.array([4/9, 1/9, 1/9, 1/9, 1/9, 1/36, 1/36, 1/36, 1/36])

    def streaming(self, u):
        u[1,:, 1:] = u[1,:, :-1]  # vers la droite
        u[2,:-1, :] = u[2, 1:, :]  # vers le haut
        u[3,:, :-1] = u[3,:, 1:]  # vers la gauche
        u[4,1:, :] = u[4,:-1, :]  # vers le bas
        u[5,:-1, 1:] = u[5,1:, :-1]  # diagonale haut droite
        u[6,:-1, :-1] = u[6,1:, 1:]  # diagonale haut gauche
        u[7,1:, :-1] = u[7,:-1, 1:]  # diagonale bas gauche
        u[8,1:, 1:] = u[8,:-1, :-1]  # diagonale bas droite
        return u
     
    @abstractmethod
    def conditions_bords(self): # méthode à définir dans les classes filles (les 2 schémas)
        pass

    def calcul_vitesse(self, u, rho):
        v = np.sum(list(map(lambda k: np.einsum('ij,k->kij', u[k],self.c[k]), range(1,9))), axis = 0)
        return v/rho

class LBM_F(LBM): # classe LBM_F pour le schéma 1 hérite de LBM
    def __init__(self, m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D):
        super().__init__(m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D)
        self.f = np.ones((9, m, m)) 
        self.init_f()
    
    def init_f(self):
        self.f[0, :, : ] = 4/9*self.rho
        self.f[1:5, :, :] = (1/9) * self.rho
        self.f[5:, :, :] = (1/36) * self.rho

    def calcul_feq(self, v, rho):
        vx, vy = v[0], v[1]
        # Calcul des fonctions d'équilibre
        #feq = np.array(list(map(lambda k: self.t[k] * rho * (1 + (self.c[k, 0] * vx + self.c[k, 1] * vy) / self.c_s2 +
        #                                             0.5 * ((self.c[k, 0] * vx + self.c[k, 1] * vy) ** 2) / self.c_s2 ** 2 -
        #     
        #                                         0.5 * (vx ** 2 + vy ** 2) / self.c_s2), range(9))))
        feq = np.zeros((9, m, m))
        for k in range(9):
            cu =  self.c[k,0] * vx + self.c[k,1] * vy
            feq[k] = self.t[k] * rho * (1 + cu / self.c_s2 + 0.5 * (cu**2) / (self.c_s2)**2 - 0.5 * (vx**2 + vy**2) / (self.c_s2))
        return feq
    
    def conditions_bords(self):
        # bord gauche
        self.f[1, :, 0] = self.f[3, :, 0]
        self.f[5, :, 0] = self.f[7, :, 0]
        self.f[8, :, 0] = self.f[6, :, 0]
        # bord droit
        self.f[3, :, -1] = self.f[1, :, -1]
        self.f[7, :, -1] = self.f[5, :, -1]
        self.f[6, :, -1] = self.f[8, :, -1]
        # bord inférieur
        self.f[2, -1, :] = self.f[4, -1, :]
        self.f[5, -1, :] = self.f[7, -1, :]
        self.f[6, -1, :] = self.f[8, -1, :]

        # bord supérieur (conditions de Zhu et He)
        rho_top = self.f[0, 0, :] + self.f[1, 0, :] + self.f[3, 0, :] + 2 * (self.f[2, 0, :] + self.f[5, 0, :] + self.f[6, 0, :])
        self.f[4, 0, :] = self.f[2, 0, :]  
        self.f[7, 0, :] = self.f[5, 0, :] - (1/6) * rho_top * self.v0 
        self.f[8, 0, :] = self.f[6, 0, :] + (1/6) * rho_top * self.v0  
        self.v[0, 0, :] = self.v0  
        self.v[1,0, :] = 0   
        return self.f

class LBM_G(LBM): # Classe pour le schéma 2
    def __init__(self, m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D):
        super().__init__(m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D)
        self.g = np.zeros((9, m, m))  


    def calcul_geq(self, v, chi):
        vx, vy = v[0], v[1]

        # calcul des fonctions d'équilibre pour le schéma 2
        # geq = np.array(list(map(lambda k: self.t[k] * chi * (1 + ((vx * self.c[k, 0] + vy * self.c[k, 1]) / self.c_s2)), range(9))))
        geq = np.zeros((9, m, m))
        for k in range(9):
            geq[k] = self.t[k] * chi * (1 + ((vx * self.c[k][0] + vy * self.c[k][1]) / (self.c_s2)))
        return geq

    def conditions_bords(self, chi):
        # bord gauche
        self.g[1, :, 0] = -self.g[3, :, 0]
        self.g[5, :, 0] = -self.g[7, :, 0]
        self.g[8, :, 0] = -self.g[6, :, 0]
        chi[: , 0] = 0
        # bord droit
        self.g[6, :, -1] = -self.g[8, :, -1]
        self.g[3, :, -1] = -self.g[1, :, -1]
        self.g[7, :, -1] = -self.g[5, :, -1]
        chi[: , -1] = 0
        # bord supérieur
        self.g[8, 0, :] = self.U0 / 36 + self.U0 / 36 - self.g[6, 0, :]
        self.g[7, 0, :] = self.U0 / 36 + self.U0 / 36 - self.g[5, 0, :]
        self.g[4, 0, :] = self.U0 / 9 + self.U0 / 9 - self.g[2, 0, :]
        chi[0, :] = self.U0 
        # bord inférieur (condition de gradient nul)

        for k in range(9):
            self.g[k, -1, :] = self.g[k, -2, :]
        return self.g, chi

class Simulation(QObject):

    def __init__(self, iterations, m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D):
        # Initialisation des paramètres de la simulation
        self.iterations = iterations
        # Initialisation des modèles LBM
        self.lbmf = LBM_F(m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D)
        self.lbmg = LBM_G(m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D)
        # Initialisation des données pour l'animation
        self.CHI = []
        self.chi = self.lbmf.chi
        self.v = self.lbmf.v
        self.rho = self.lbmf.rho 

        self.vitesse_animation = 1
        self.V = []
        self.RHO = []

        self.intG = None
    def set_intG(self, intG):
        self.intG = intG

    def set_nb_iteration(self, steps):
        self.iterations = steps

    def set_vitesse(self, v):
        self.vitesse_animation = v

    def simul(self): # Boucle principale de simulation
        for i in range(self.iterations):
            self.intG.progress_bar.setValue(i+1)
            #↕self.progression_changed.emit()
            # Données de collision pour les fk
            feq = self.lbmf.calcul_feq(self.v, self.rho)
            # Collision pour les fk
            for k in range(9):
                self.lbmf.f[k] = self.lbmf.f[k] * (1 - self.lbmf.omega_m) + self.lbmf.omega_m * feq[k]
            # Streaming pour les fk
            self.lbmf.f = self.lbmf.streaming(self.lbmf.f)
            # Conditions aux bords pour les fk
            self.lbmf.f = self.lbmf.conditions_bords()
            # Calcul de ρ et de v
            self.rho = np.sum(self.lbmf.f, axis=0)
            self.v = self.lbmf.calcul_vitesse(self.lbmf.f, self.rho)
            # Données de collision pour les gk
            geq= self.lbmg.calcul_geq(self.v, self.chi)
            # Collision pour les gk
            for k in range(9):
                self.lbmg.g[k] = self.lbmg.g[k] * (1 - self.lbmg.omega_s) + self.lbmg.omega_s * geq[k]
            # Streaming pour les gk
            self.lbmg.g = self.lbmg.streaming(self.lbmg.g)
            # Conditions aux bords pour les gk
            self.lbmg.g, self.chi = self.lbmg.conditions_bords(self.chi)
            # Calcul de χ 
            self.chi = np.sum(self.lbmg.g, axis=0)

            # Stockage des résultats pour l'animation
            if i %self.vitesse_animation == 0: # 

                self.CHI.append(self.chi)

                self.V.append(self.v)

                self.RHO.append(self.rho)

        return self.CHI, self.V, self.RHO
    

if __name__ == "__main__":
    m = 100
    delta = 1
    tau = 1
    rho0 = 5
    chi0 = 0
    v0 = 0.2
    U0 = 1
    nu = 0.02
    Re = 1000
    Pr = 0.71
    D = 0.02817

    app = QApplication(sys.argv)
    interface = InterfaceGraphique(m, delta, tau, rho0, chi0, v0, U0, nu, Re, Pr, D)
    interface.setGeometry(100, 100, 400, 400)
    interface.show()
    sys.exit(app.exec_())
