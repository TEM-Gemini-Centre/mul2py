% Function for calculating the Scherzer defocus in [Å] given an acceleration voltage in [kV] and third order spherical aberration Cs in [mm].
function [defocus] = scherzer_defocus(E0, Cs)
    Cs = Cs * 1E-3*1E10; %Convert from mm to Å
    defocus = -sqrt(Cs * energy2wavelength(E0) * 4/3);
end