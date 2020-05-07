% Function for converting acceleration voltage in [kV] of an electron to electron wavelength in [Å]
function [wavelength] = energy2wavelength(E0)
    e=1.602E-19;    %Electron charge [C]
    m0=9.109E-31;   %Rest mass of electron [kg]
    h=6.626E-34;    %Plancks constant [N m s]
    c=2.998E8;      %Speed of light in vacuum [m]

    E0 = E0*1E3;    %kV to V

    wavelength =  h / (sqrt(2 * m0 * e * E0 * (1. + e * E0 / (2. * m0 * c^2)))); %wavelength in m
    wavelength = wavelength * 1E10; %m to Å
end