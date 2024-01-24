import numpy as np
import statsmodels.api as sm
import miepython
import matplotlib.pyplot as plt
import math

sca_450_pm_1 = 11.3
sca_550_pm_1 = 7.1
sca_700_pm_1 = 3.8

sca_450_pm_10 = 14.2
sca_550_pm_10 = 9.8
sca_700_pm_10 = 6.3

bsca_450_pm_1 = 1.4
bsca_550_pm_1 = 1.1
bsca_700_pm_1 = 0.8

bsca_450_pm_10 = 1.8
bsca_550_pm_10 = 1.4
bsca_700_pm_10 = 1.2

abs_370_pm_1 = 1.6
abs_520_pm_1 = 1.2
abs_880_pm_1 = 0.7

abs_370_pm_10 = 1.9
abs_520_pm_10 = 1.4
abs_880_pm_10 = 0.9

bscaf_450_pm_1 = bsca_450_pm_1/sca_450_pm_1
bscaf_550_pm_1 = bsca_550_pm_1/sca_550_pm_1
bscaf_700_pm_1 = bsca_700_pm_1/sca_700_pm_1

bscaf_450_pm_1 = round(bscaf_450_pm_1, 2)
bscaf_550_pm_1 = round(bscaf_550_pm_1, 2)
bscaf_700_pm_1 = round(bscaf_700_pm_1, 2)

print('bscaf_450_pm_1:', bscaf_450_pm_1)
print('bscaf_550_pm_1:', bscaf_550_pm_1)
print('bscaf_700_pm_1:', bscaf_700_pm_1)
print()

bscaf_450_pm_10 = bsca_450_pm_10/sca_450_pm_10
bscaf_550_pm_10 = bsca_550_pm_10/sca_550_pm_10
bscaf_700_pm_10 = bsca_700_pm_10/sca_700_pm_10

bscaf_450_pm_10 = round(bscaf_450_pm_10, 2)
bscaf_550_pm_10 = round(bscaf_550_pm_10, 2)
bscaf_700_pm_10 = round(bscaf_700_pm_10, 2)

print('bscaf_450_pm_10:', bscaf_450_pm_10)
print('bscaf_550_pm_10:', bscaf_550_pm_10)
print('bscaf_700_pm_10:', bscaf_700_pm_10) 
print() 

absorption_values = np.array([1.6, 1.2, 0.7])
wavelengths = np.array([370, 520, 880])

ln_absorption_values = np.log(absorption_values)
ln_wavelengths = np.log(wavelengths)

ln_wavelengths_with_const = sm.add_constant(ln_wavelengths)

model = sm.OLS(ln_absorption_values, ln_wavelengths_with_const)
results = model.fit()
aae = -results.params[1]

print('The absorption Ångström exponent is:', aae)
print()

abs_550_pm_1 = abs_520_pm_1/ (np.exp(np.log(520/550)*(-aae)))
abs_550_pm_1 = round(abs_550_pm_1, 2)

print('abs_550_pm_1:', abs_550_pm_1)

absorption_values = np.array([1.9, 1.4, 0.9])
wavelengths = np.array([370, 520, 880])

ln_absorption_values = np.log(absorption_values)
ln_wavelengths = np.log(wavelengths)

ln_wavelengths_with_const = sm.add_constant(ln_wavelengths)

model = sm.OLS(ln_absorption_values, ln_wavelengths_with_const)
results = model.fit()
aae = -results.params[1]

print('The absorption Ångström exponent is:', aae)
print()

abs_550_pm_10 = abs_520_pm_10/ (np.exp(np.log(520/550)*(-aae)))
abs_550_pm_10 = round(abs_550_pm_10, 2)

print('abs_550_pm_10:', abs_550_pm_10)
print()

ext_550_pm_1 = abs_550_pm_1 + sca_550_pm_1
print('ext_550_pm_1:', ext_550_pm_1)
print()

ext_550_pm_10 = abs_550_pm_10 + sca_550_pm_10
print('ext_550_pm_10:', ext_550_pm_10)
print()

ssa_550_pm_1 = sca_550_pm_1/ext_550_pm_1
ssa_550_pm_1 = round(sca_550_pm_1, 2)
print('ssa_550_pm_1:', ssa_550_pm_1)
print()

ssa_550_pm_10 = sca_550_pm_10/ext_550_pm_10
ssa_550_pm_10 = round(sca_550_pm_10, 2)
print('ssa_550_pm_10:', ssa_550_pm_10)
print()

plt.figure(figsize = (30, 15))

wavelength = 550
min_diameter = 10
max_diameter = 10000

def calculate_size_parameter(radius, wavelength):
    return 2 * np.pi * radius / wavelength

m_water = 1.333
N_water = 600
radius_water = np.linspace(min_diameter / 2, max_diameter / 2, N_water)
x_water = calculate_size_parameter(radius_water, wavelength)
qext_water, qsca_water, qback_water, g_water = miepython.mie(m_water, x_water)

m_sulphuric_acid = 1.431
N_sulphuric_acid = 600
radius_sulphuric_acid = np.linspace(min_diameter / 2, max_diameter / 2, N_sulphuric_acid)
x_sulphuric_acid = calculate_size_parameter(radius_sulphuric_acid, wavelength)
qext_sulphuric_acid, qsca_sulphuric_acid, qback_sulphuric_acid, g_sulphuric_acid = miepython.mie(m_sulphuric_acid, x_sulphuric_acid)

m_ammonium_sulphate = 1.521
N_ammonium_sulphate = 600
radius_ammonium_sulphate = np.linspace(min_diameter / 2, max_diameter / 2, N_ammonium_sulphate)
x_ammonium_sulphate = calculate_size_parameter(radius_ammonium_sulphate, wavelength)
qext_ammonium_sulphate, qsca_ammonium_sulphate, qback_ammonium_sulphate, g_ammonium_sulphate = miepython.mie(m_ammonium_sulphate, x_ammonium_sulphate)

plt.plot(x_water, qsca_water, label = 'H$_2$O')
plt.plot(x_sulphuric_acid, qsca_sulphuric_acid, label = 'H$_2$SO$_4$')
plt.plot(x_ammonium_sulphate, qsca_ammonium_sulphate, label = '(NH$_4$)$_2$SO$_4$')

plt.xscale('log')
plt.xlabel('x')
plt.ylabel('Q$_s$')
plt.title('Scattering efficiencies at 550 nm')
plt.legend()
plt.show()
print()

plt.figure(figsize = (30, 15))

particle_radius = 200/2 
wavelengths = np.arange(300, 700, 10)

n_ammonium_sulphate = 1.521

n_ammonium_sulphate_list = [0.001, 0.01, 0.1, 0.5]

qs_list = []
qa_list = []
ssa_list = []

for n in n_ammonium_sulphate_list:
    qs_n_list = []
    qa_n_list = []

    for wavelength in wavelengths:
        x = (2 * np.pi * particle_radius)/wavelength

        qext, qsca, qback, g = miepython.mie(n_ammonium_sulphate - 1j * n, x)
        qs_n_list.append(qsca)
        qa_n_list.append(qext - qsca)

    qs_list.append(qs_n_list)
    qa_list.append(qa_n_list)

    ssa_list.append(np.array(qs_n_list)/np.array(qa_n_list))

qs_values = np.array(qs_list)
qa_values = np.array(qa_list)
ssa_values = np.array(ssa_list)

for i, n_ammonium_sulphate in enumerate(n_ammonium_sulphate_list):
    plt.plot(wavelengths, qs_values[i], label = f'n(NH$_4$)$_2$SO$_4$ = {n_ammonium_sulphate}')

plt.xlabel('Wavelength (nm)')
plt.ylabel('$Q_{\\mathrm{scat}}$')
plt.title('Scattering efficiencies of 200 nm (NH$_4$)$_2$SO$_4$ particles')
plt.legend()
plt.show()
print()

plt.figure(figsize = (30, 15))

for i, n_ammonium_sulphate in enumerate(n_ammonium_sulphate_list):
    plt.plot(wavelengths, ssa_list[i], label = f'n(NH$_4$)$_2$SO$_4$ = {n_ammonium_sulphate}')
    
legend_label = ', '.join([f'n(NH$_4$)$_2$SO$_4$ = {n}' for n in n_ammonium_sulphate_list])

plt.xlabel('Wavelength (nm)')
plt.ylabel('$Q_{\\mathrm{scat}}$')
plt.title('Scattering efficiencies of 200 nm (NH$_4$)$_2$SO$_4$ particles')
plt.legend()
plt.show()
print()

sigma_ext_ambair = 13.2 * 10**(-6)

visibility = -math.log(0.02) / sigma_ext_ambair

print('Visibility:', visibility, 'm')
print()

visibility = 5000
particle_diameter = 100 * 10**(-9)

N = -math.log(0.02) / (visibility * math.pi * (particle_diameter/2)**2)

print(f'Particle Number Concentration (N): {N:.2e} particles/m^3')
print()

visibility_km = 80
visibility_m = visibility_km * 1000

sigma_ext_total = -math.log(0.02) / visibility_m

print(f"Extinction Coefficient (sigma_ext_total): {sigma_ext_total:.2e} m^-1")


















