class strain:
    def __init__(self, h, t, f22, df22dt, dt, i1, i2):
        self.h = h
        self.t = t
        self.f22 = f22
        self.df22dt = df22dt
        self.dt = dt
        self.i1 = i1
        self.i2 = i2

    def initialize(self, G, c, M, r):
        h.time = h.time * G * (M/(c**3))
        h = h * (M/r) * (G/(c**2))
        t = h.t
        h22 = h.data[:, h.index(2,2)]
        A22 = np.abs(h22)
        phi22 = np.unwrap(np.angle(h22))
        phi22_t = scipy.interpolate.CubicSpline(t, phi22)
        f22 = phi22_t.derivative()(t) / (2*np.pi)
        df22dt = phi22_t.derivative(2)(t) / (2*np.pi)
        dt = 1e-4 #np.min(np.diff(h.t))
        i1 = h.index_closest_to(0.5)
        i2 = h.max_norm_index()