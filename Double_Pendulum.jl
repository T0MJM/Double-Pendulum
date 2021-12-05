using Plots
using DifferentialEquations
using Dierckx

pyplot()

m1 = 1;
m2 = 2;
L1 = 1;
l2 = 2;
g = 9.81;

u0 = [pi/1.6; 0; pi/1.8; 0]; #Initial conditions

tfinal = 100.0               #Final time

#Pendulum ODE
function dbl_pendulum(du,u,p,t)

    m1 = p[1];
    m2 = p[2];
    L1 = p[3];
    L2 = p[4];
    g = p[5];

    c = cos(u[1]-u[3])
    s = sin(u[1]-u[3])

    du[1] = u[2]
    du[2] = (m2*g*sin(u[3])*c - m2*s*(L1*c*u[2]^2 + L2*u[4]^2) - (m1+m2)*g*sin(u[1])) / (L1*(m1+m2*s^2))
    du[3] = u[4];
    du[4] = ((m1+m2)*(L1*u[2]^2*s - g*sin(u[3]) + g*sin(u[1])*c) + m2*L2*u[4]^2*s*c) / (L2*(m1 + m2*s^2))

end

p = [m1;  m2;  L1;  L2;  g];
tspan = (0.0,tfinal);  # Time span (limits). The actual time variable is automatically set by solve().
prob = ODEProblem(double_pendulum,u0,tspan,p);
sol = solve(prob,Vern7(),reltol=1e-6);

#Extract variables

tm = sol.t;

# Mapping from polar to Cartesian

x1 = L1*sin.(sol[1,:]);          # First Pendulum
y1 = -L1*cos.(sol[1,:]);

x2 = x1 + L2*sin.(sol[3,:]);     # Second Pendulum
y2 = y1 - L2*cos.(sol[3,:]);


#Interpolation

dt = 0.05;
t_u = 0:dt:tfinal;      # uniformly spaced time variable

sp_x1 = Spline1D(tm, x1);
sp_y1 = Spline1D(tm, y1);
sp_x2 = Spline1D(tm, x2);
sp_y2 = Spline1D(tm, y2);


# Interpolated variables
x1_u = sp_x1(t_u);
y1_u = sp_y1(t_u);
x2_u = sp_x2(t_u);
y2_u = sp_y2(t_u);