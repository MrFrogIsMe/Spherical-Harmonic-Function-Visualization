from vpython import *
import sympy as sym

# x_sym is the symbol for the associated Legendre polynomial
x_sym = sym.symbols('x')

# define the associated Legendre polynomial
def P_n(n):
    return 1/(2**n * sym.factorial(n)) * sym.diff((x_sym**2 - 1)**n, x_sym, n)

# define the associated Legendre polynomial derivative
def P_nm(n, m):
    return (1 - x_sym**2)**(n/2) * sym.diff(P_n(n), x_sym, sym.Abs(m))

# define the spherical harmonic function
def Y_nm(n, m, theta, phi):
    pre_factor = (-1)**m * sym.sqrt((2*n+1)/(4*sym.pi)*sym.factorial(n-abs(m))/sym.factorial(n+abs(m)))
    theta_factor = P_nm(n, m).subs(x_sym, sym.cos(theta))
    phase_factor = sym.exp(1j*m*phi)
    return pre_factor * theta_factor * phase_factor

title = '''Spherical Harmonic Function

Rotate the camera view: drag with the right mouse button (or Ctrl-drag left button).
Zoom: drag with left and right mouse buttons (or Alt/Option-drag or scroll wheel).
Pan: Shift-drag.
'''
scene = canvas(title=title, width=800, height=600, x=0, y=0, center=vec(0, 0, 0), background=vec(0, 0.6, 0.6))

# set z axis to be upward
scene.camera.rotate(angle=pi/2, axis=vec(1, 0, 0), origin=vec(0, 0, 0))
scene.camera.rotate(angle=pi/2, axis=vec(0, 0, 1), origin=vec(0, 0, 0))
scene.center = vec(0, 0, 0)

# draw axis for x, y, z
axis_x = arrow(pos=vec(0, 0, 0), axis=vec(1, 0, 0)*3, color=vec(1, 0, 0))
axis_y = arrow(pos=vec(0, 0, 0), axis=vec(0, 1, 0)*3, color=vec(0, 1, 0))
axis_z = arrow(pos=vec(0, 0, 0), axis=vec(0, 0, 1)*3, color=vec(0, 0, 1))

# set the width of the axis
axis_x.shaftwidth = 0.02
axis_y.shaftwidth = 0.02
axis_z.shaftwidth = 0.02

# set the label of the axis
axis_x_label = label(pos=axis_x.pos+axis_x.axis, text='x', xoffset=20, yoffset=20, space=30, height=16, border=4, font='sans')
axis_y_label = label(pos=axis_y.pos+axis_y.axis, text='y', xoffset=20, yoffset=20, space=30, height=16, border=4, font='sans')
axis_z_label = label(pos=axis_z.pos+axis_z.axis, text='z', xoffset=20, yoffset=20, space=30, height=16, border=4, font='sans')

# set up the variables for the spherical harmonic function
n = 1
m = 1
r = 1
info_label = label(pos=vec(0, 0, 0), text=f'n = {n}, m = {m}', xoffset=20, yoffset=20, space=30, height=16, border=4, font='sans')
step = 0.06
planes = [[box(pos=vec(0, 0, 0), size=vec(0.07, 0.001, 0.07), color=vec(1, 1, 1)) for theta in arange(0, pi, step)] for phi in arange(0, 2*pi, step)]
for i, phi in enumerate(arange(0, 2*pi, step)):
    for j, theta in enumerate(arange(0, pi, step)):
        x = sin(theta) * cos(phi)
        y = sin(theta) * sin(phi)
        z = cos(theta)
        planes[i][j].pos = vec(x, y, z) * r
        planes[i][j].up = vec(x, y, z)
        
while True:
    rate(1000)
    info_label.text = 'l = %d, m = %d' % (n, m)
    f = Y_nm(n, m, sym.symbols('theta'), sym.symbols('phi'))
    for i, phi in enumerate(arange(0, 2*pi, step)):
        for j, theta in enumerate(arange(0, pi, step)):
            # solve the spherical harmonic function
            solve = f.subs(sym.symbols('theta'), theta).subs(sym.symbols('phi'), phi).evalf()
            solve = sym.re(solve)

            # convert the spherical coordinate to Cartesian coordinate
            x = sin(theta) * cos(phi)
            y = sin(theta) * sin(phi)
            z = cos(theta)

            # set the color of the surface
            # if the real part of the spherical harmonic function is positive, set the color to white, otherwise set the color to black
            planes[i][j].color = vec(1, 1, 1) if solve > 0 else vec(0, 0, 0)

    # wait for the user to enter the value of l and m
    n, m = map(int, input('Enter n, m: ').split())
