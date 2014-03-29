function [] = fdtd_mat2h5( input_file, h5_file )
% NOTE: This program is free software; you can redistribute it
% and/or modify it under the terms of the GNU General Public 
% License as published by the Free Software Foundation; either
% version 3, or (at you option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.

if exist(input_file, 'file') ~= 2;
    error(['file "' input_file '" does not exsit.']); % remove the file to avoid overwrite
end

if exist(h5_file, 'file') == 2;
    delete(h5_file); % remove the file to avoid overwrite
end
load(input_file);

% create groups
h5create(h5_file, '/settings/dummy', 1);
h5create(h5_file, '/excitations/dummy', 1);
h5create(h5_file, '/materials/dummy', 1);
h5create(h5_file, '/boundaries/dummy', 1);
h5create(h5_file, '/outputs/dummy', 1);
h5create(h5_file, '/DFT/dummy', 1);

%=================================
% settings
%=================================
if exist('verbose_level', 'var') == 0; 
    verbose_level = 1;
end
if exist('simulation_mode', 'var') == 0; 
    error([input_file ': scalar "simulation_mode" does not exsit']); 
end
if exist('NUMBER_OF_TIMESTEPS', 'var') == 0; 
    error([input_file ': scalar "NUMBER_OF_TIMESTEPS" does not exsit']); 
end
if exist('domain_size_x', 'var') == 0; 
    error([input_file ': scalar "domain_size_x" does not exsit']); 
end
if exist('domain_size_y', 'var') == 0; 
    error([input_file ': scalar "domain_size_y" does not exsit']); 
end
if exist('domain_size_z', 'var') == 0; 
    error([input_file ': scalar "domain_size_x" does not exsit']); 
end
h5writeatt(h5_file,'/settings', 'starting_index', int32(1));
h5writeatt(h5_file,'/settings', 'verbosity', int32(verbose_level));
h5writeatt(h5_file,'/settings', 'mode', int32(simulation_mode));
h5writeatt(h5_file,'/settings', 'number_of_timesteps', int32(NUMBER_OF_TIMESTEPS));
h5writeatt(h5_file,'/settings', 'domain_size_x', int32(domain_size_x));
h5writeatt(h5_file,'/settings', 'domain_size_y', int32(domain_size_y));
h5writeatt(h5_file,'/settings', 'domain_size_z', int32(domain_size_z));

if exist('d_t', 'var') == 0; 
    error([input_file ': scalar "d_t" does not exsit']); 
end
if exist('d_x', 'var') == 0; 
    error([input_file ': scalar "d_x" does not exsit']); 
end
if exist('d_y', 'var') == 0; 
    error([input_file ': scalar "d_y" does not exsit']); 
end
if exist('d_z', 'var') == 0; 
    error([input_file ': scalar "d_z" does not exsit']); 
end
h5writeatt(h5_file,'/settings', 'd_t', d_t);
h5writeatt(h5_file,'/settings', 'd_x', d_x);
h5writeatt(h5_file,'/settings', 'd_y', d_y);
h5writeatt(h5_file,'/settings', 'd_z', d_z);

if exist('C0', 'var') == 0; 
    C0 = 299792458.0;
end
if exist('MU0', 'var') == 0; 
    MU0 = 4.0 * pi * 1.0e-7;
end
if exist('EPSILON0', 'var') == 0; 
    EPSILON0 = 1.0 / (C0^2 * MU0);
end
h5writeatt(h5_file,'/settings', 'C0', C0);
h5writeatt(h5_file,'/settings', 'EPSILON0', EPSILON0);
h5writeatt(h5_file,'/settings', 'MU0', MU0);


%=================================
% dft
%=================================
if exist('ENABLE_DFT', 'var') == 0; 
    ENABLE_DFT = 0;
end
h5writeatt(h5_file,'/DFT', 'enable_DFT', int32(ENABLE_DFT));
if int32(ENABLE_DFT) == 1;
    if exist('NUMBER_OF_FREQUENCIES', 'var') == 0; 
        error([h5_file ': scalar "NUMBER_OF_FREQUENCIES" does not exsit']); 
    end
    if exist('frequencies', 'var') == 0; 
        error([h5_file ': vector "frequencies" does not exsit']); 
    end 
    h5writeatt(h5_file,'/DFT', 'number_of_frequencies', int32(NUMBER_OF_FREQUENCIES));
    h5create(h5_file, '/DFT/frequencies', [NUMBER_OF_FREQUENCIES]);
    h5write(h5_file, '/DFT/frequencies', frequencies);
end

%=================================
% input ports
%=================================
if exist('NUMBER_OF_INPUTPORTS', 'var') == 0; 
    NUMBER_OF_INPUTPORTS = 0;
end
h5writeatt(h5_file,'/excitations', 'number_of_ports', int32(NUMBER_OF_INPUTPORTS));
for port_index = 1:NUMBER_OF_INPUTPORTS;
    var_name = sprintf('inputport_%02d_number_of_points', port_index);
    attr_name = sprintf('port_%02d_number_of_points', port_index);
    var_dim = eval(var_name);
    if exist(var_name, 'var') == 0; 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    end
    h5writeatt(h5_file,'/excitations', attr_name, int32(eval(var_name)));

    var_name = sprintf('inputport_%02d_polarization', port_index);
    attr_name = sprintf('port_%02d_polarization', port_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    end
    h5writeatt(h5_file,'/excitations', attr_name, int32(eval(var_name))); 

    var_name = sprintf('inputport_%02d_points', port_index);
    dset_name = sprintf('/excitations/port_%02d_points', port_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [3 var_dim], 'Datatype','int32');
    h5write (h5_file, dset_name, int32(eval(var_name))');

    var_name = sprintf('inputport_%02d_distribution', port_index);
    dset_name = sprintf('/excitations/port_%02d_distribution', port_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [var_dim]);
    h5write(h5_file, dset_name, eval(var_name));

    var_name = sprintf('inputport_%02d_signal', port_index);
    dset_name = sprintf('/excitations/port_%02d_signal', port_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [NUMBER_OF_TIMESTEPS]);
    h5write(h5_file, dset_name, eval(var_name));
end

%=================================
% point sources
%=================================
if exist('NUMBER_OF_POINT_SOURCES', 'var') == 0; 
    NUMBER_OF_POINT_SOURCES = 0;
end
h5writeatt(h5_file,'/excitations', 'number_of_point_sources', int32(NUMBER_OF_POINT_SOURCES));
if NUMBER_OF_POINT_SOURCES > 0;
    if exist('source_points', 'var') == 0; 
        error([input_file ': matrix ' 'source_points' ' does not exsit']); 
    end
    if exist('source_polarizations', 'var') == 0; 
        error([input_file ': matrix ' 'source_polarizations' ' does not exsit']); 
    end
    if exist('source_signals', 'var') == 0; 
        error([input_file ': matrix ' 'source_signals' ' does not exsit']); 
    end
    h5create(h5_file, '/excitations/source_points', [3 NUMBER_OF_POINT_SOURCES], 'Datatype','int32');
    h5create(h5_file, '/excitations/source_polarizations', [NUMBER_OF_POINT_SOURCES], 'Datatype','int32');
    h5create(h5_file, '/excitations/source_signals', [NUMBER_OF_POINT_SOURCES NUMBER_OF_TIMESTEPS]);
    h5write (h5_file, '/excitations/source_points', int32(source_points'));
    h5write (h5_file, '/excitations/source_polarizations', int32(source_polarizations));
    h5write (h5_file, '/excitations/source_signals', source_signals');
end

%=================================
% hard sources
%=================================
if exist('NUMBER_OF_HARDS', 'var') == 0; 
    NUMBER_OF_HARDS = 0;
end
h5writeatt(h5_file,'/excitations', 'number_of_hards', int32(NUMBER_OF_HARDS));
if NUMBER_OF_HARDS > 0;
    if exist('hard_points', 'var') == 0; 
        error([input_file ': matrix ' 'hard_points' ' does not exsit']); 
    end
    if exist('hard_polarizations', 'var') == 0; 
        error([input_file ': matrix ' 'hard_polarizations' ' does not exsit']); 
    end
    if exist('hard_signals', 'var') == 0; 
        error([input_file ': matrix ' 'hard_signals' ' does not exsit']); 
    end
    h5create(h5_file, '/excitations/hard_points', [3 NUMBER_OF_HARDS], 'Datatype','int32');
    h5create(h5_file, '/excitations/hard_polarizations', [NUMBER_OF_HARDS], 'Datatype','int32');
    h5create(h5_file, '/excitations/hard_signals', [NUMBER_OF_HARDS NUMBER_OF_TIMESTEPS]);
    h5write (h5_file, '/excitations/hard_points', int32(hard_points'));
    h5write (h5_file, '/excitations/hard_polarizations', int32(hard_polarizations));
    h5write (h5_file, '/excitations/hard_signals', hard_signals');
end

%=================================
% animation planes
%=================================
if exist('NUMBER_OF_ANIMATION_PLANES', 'var') == 0; 
    NUMBER_OF_ANIMATION_PLANES = 0;
end
h5writeatt(h5_file,'/outputs', 'number_of_field_planes', int32(NUMBER_OF_ANIMATION_PLANES));
for plane_index = 1:NUMBER_OF_ANIMATION_PLANES;
    var_name  = sprintf('plane_%02d_polarization', plane_index);
    attr_name = sprintf('plane_%02d_polarization', plane_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    end
    h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));

    var_name  = sprintf('plane_%02d_dim', plane_index);
    attr_name = sprintf('plane_%02d_dim', plane_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    end
    h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));

    var_name  = sprintf('plane_%02d_slice_index', plane_index);
    attr_name = sprintf('plane_%02d_slice_index', plane_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    end
    h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));

    var_name  = sprintf('plane_%02d_x_start', plane_index);
    attr_name = sprintf('plane_%02d_x_start', plane_index);
    if exist(var_name, 'var') == 0 & eval([var_name(1:8) '_dim ~= 1']); 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    else
        if exist(var_name, 'var') == 1;
            h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));
        end
    end

    var_name  = sprintf('plane_%02d_y_start', plane_index);
    attr_name = sprintf('plane_%02d_y_start', plane_index);
    if exist(var_name, 'var') == 0 & eval([var_name(1:8) '_dim ~= 2']); 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    else
        if exist(var_name, 'var') == 1;
            h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));
        end
    end

    var_name  = sprintf('plane_%02d_z_start', plane_index);
    attr_name = sprintf('plane_%02d_z_start', plane_index);
    if exist(var_name, 'var') == 0 & eval([var_name(1:8) '_dim ~= 3']); 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    else
        if exist(var_name, 'var') == 1;
            h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));
        end
    end

    var_name  = sprintf('plane_%02d_length_x', plane_index);
    attr_name = sprintf('plane_%02d_length_x', plane_index);
    if exist(var_name, 'var') == 0 & eval([var_name(1:8) '_dim ~= 1']); 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    else
        if exist(var_name, 'var') == 1;
            h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));
        end
    end

    var_name  = sprintf('plane_%02d_length_y', plane_index);
    attr_name = sprintf('plane_%02d_length_y', plane_index);
    if exist(var_name, 'var') == 0 & eval([var_name(1:8) '_dim ~= 2']); 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    else
        if exist(var_name, 'var') == 1;
            h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));
        end
    end

    var_name  = sprintf('plane_%02d_length_z', plane_index);
    attr_name = sprintf('plane_%02d_length_z', plane_index);
    if exist(var_name, 'var') == 0 & eval([var_name(1:8) '_dim ~= 3']); 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    else
        if exist(var_name, 'var') == 1;
            h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));
        end
    end
end

if exist('NUMBER_OF_OUTPUTPORTS', 'var') == 0; 
    NUMBER_OF_OUTPUTPORTS = 0;
end
for port_index = 1:NUMBER_OF_OUTPUTPORTS;

end


%=================================
% output ports
%=================================
if exist('NUMBER_OF_OUTPUTPORTS', 'var') == 0; 
    NUMBER_OF_OUTPUTPORTS = 0;
end
h5writeatt(h5_file,'/outputs', 'number_of_ports', int32(NUMBER_OF_OUTPUTPORTS));
for port_index = 1:NUMBER_OF_OUTPUTPORTS;
    var_name = sprintf('outputport_%02d_write_timedomain', port_index);
    attr_name = sprintf('port_%02d_write_timedomain', port_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    end
    h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));

    var_name = sprintf('outputport_%02d_number_of_points', port_index);
    attr_name = sprintf('port_%02d_number_of_points', port_index);
    var_dim = eval(var_name);
    if exist(var_name, 'var') == 0; 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    end
    h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));

    var_name = sprintf('outputport_%02d_polarization', port_index);
    attr_name = sprintf('port_%02d_polarization', port_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': scalar ' var_name ' does not exsit']); 
    end
    h5writeatt(h5_file,'/outputs', attr_name, int32(eval(var_name)));

    dset_name = sprintf('/outputs/port_%02d_points', port_index);
    var_name = sprintf('outputport_%02d_points', port_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [3 var_dim], 'Datatype','int32');
    h5write (h5_file, dset_name, int32(eval(var_name))');
end

%=================================
% boundaries
%=================================
if exist('pml_type', 'var') == 0; 
    pml_type = 0;
end
if exist('abc_thickness', 'var') == 0; 
    abc_thickness = 200;
end
h5writeatt(h5_file,'/boundaries', 'abc_thickness', int32(abc_size));
if exist('pml_grading_order', 'var') == 0; 
    pml_grading_order = 3.0;
end
if exist('pml_sigma_max', 'var') == 0; 
    pml_sigma_max = 0.8 * (pml_grading_order + 1) / (377.0 * mean([d_x d_y d_z]));
end
if exist('pml_kappa_max', 'var') == 0; 
    pml_kappa_max = 25.0;
end
if exist('pml_alpha_max', 'var') == 0; 
    pml_alpha_max = 3.0;
end
h5writeatt(h5_file,'/boundaries', 'pml_type', int32(pml_type));
h5writeatt(h5_file,'/boundaries', 'pml_grading_order', pml_grading_order);
h5writeatt(h5_file,'/boundaries', 'pml_sigma_max', pml_sigma_max);
h5writeatt(h5_file,'/boundaries', 'pml_kappa_max', pml_kappa_max);
h5writeatt(h5_file,'/boundaries', 'pml_alpha_max', pml_alpha_max);

if exist('boundary_z0', 'var') == 0; 
    error([input_file ': boundary z0 is not specified.']); 
end
if exist('boundary_z1', 'var') == 0; 
    error([input_file ': boundary z1 is not specified.']); 
end
if exist('boundary_x0', 'var') == 0; 
    error([input_file ': boundary x0 is not specified.']); 
end
if exist('boundary_x1', 'var') == 0; 
    error([input_file ': boundary x1 is not specified.']); 
end
if exist('boundary_y0', 'var') == 0; 
    error([input_file ': boundary y0 is not specified.']); 
end
if exist('boundary_y1', 'var') == 0; 
    error([input_file ': boundary y1 is not specified.']); 
end

h5writeatt(h5_file,'/boundaries', 'boundary_1', int32(boundary_z0));
h5writeatt(h5_file,'/boundaries', 'boundary_2', int32(boundary_z1));
h5writeatt(h5_file,'/boundaries', 'boundary_3', int32(boundary_x0));
h5writeatt(h5_file,'/boundaries', 'boundary_4', int32(boundary_x1));
h5writeatt(h5_file,'/boundaries', 'boundary_5', int32(boundary_y0));
h5writeatt(h5_file,'/boundaries', 'boundary_6', int32(boundary_y1));


%=================================
% material computation
%=================================
if exist('classical_mode', 'var') == 0; 
    classical_mode = 0;
end
h5writeatt(h5_file,'/materials', 'classical_mode', int32(mode_classical));

if exist('ade_sigma', 'var') == 0; 
    ade_sigma = 0;
end
h5writeatt(h5_file,'/materials', 'ade_sigma', int32(ade_sigma));

if exist('NUMBER_OF_DRUDE_POLES', 'var') == 0; 
    NUMBER_OF_DRUDE_POLES = 0;
end
if exist('NUMBER_OF_LORENTZ_POLES', 'var') == 0; 
    NUMBER_OF_LORENTZ_POLES = 0;
end
if exist('NUMBER_OF_DEBYE_POLES', 'var') == 0; 
    NUMBER_OF_DEBYE_POLES = 0;
end
h5writeatt(h5_file,'/materials', 'number_of_drude_poles', int32(NUMBER_OF_DRUDE_POLES));
h5writeatt(h5_file,'/materials', 'number_of_debye_poles', int32(NUMBER_OF_DEBYE_POLES));
h5writeatt(h5_file,'/materials', 'number_of_lorentz_poles', int32(NUMBER_OF_LORENTZ_POLES));

if exist('epsilon', 'var') == 0; 
    error([input_file ': matrix ' 'epsilon' ' does not exsit']); 
end
if exist('sigma', 'var') == 0; 
    sigma = zeros(domain_size_x, domain_size_y, domain_size_z);
end
h5create(h5_file, '/materials/epsilon', [domain_size_z domain_size_y domain_size_x]);
h5create(h5_file, '/materials/sigma', [domain_size_z domain_size_y domain_size_x]);
h5write(h5_file, '/materials/epsilon', permute(epsilon, [3,2,1]));
h5write(h5_file, '/materials/sigma', permute(sigma, [3,2,1])); 

for pole_index = 1:NUMBER_OF_DRUDE_POLES;
    dset_name = sprintf('/materials/drude_%02d_a', pole_index);
    var_name = sprintf('drude_%02d_a', pole_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [domain_size_z domain_size_y domain_size_x]);
    h5write (h5_file, dset_name, permute(eval(var_name), [3 2 1]));

    dset_name = sprintf('/materials/drude_%02d_c', pole_index);
    var_name = sprintf('drude_%02d_c', pole_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [domain_size_z domain_size_y domain_size_x]);
    h5write (h5_file, dset_name, permute(eval(var_name), [3 2 1]));
end

for pole_index = 1:NUMBER_OF_LORENTZ_POLES;
    dset_name = sprintf('/materials/lorentz_%02d_a', pole_index);
    var_name = sprintf('lorentz_%02d_a', pole_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [domain_size_z domain_size_y domain_size_x]);
    h5write (h5_file, dset_name, permute(eval(var_name), [3 2 1]));

    dset_name = sprintf('/materials/lorentz_%02d_b', pole_index);
    var_name = sprintf('lorentz_%02d_b', pole_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [domain_size_z domain_size_y domain_size_x]);
    h5write (h5_file, dset_name, permute(eval(var_name), [3 2 1]));

    dset_name = sprintf('/materials/lorentz_%02d_c', pole_index);
    var_name = sprintf('lorentz_%02d_c', pole_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [domain_size_z domain_size_y domain_size_x]);
    h5write (h5_file, dset_name, permute(eval(var_name), [3 2 1]));
end

for pole_index = 1:NUMBER_OF_DEBYE_POLES;
    dset_name = sprintf('/materials/debye_%02d_a', pole_index);
    var_name = sprintf('debye_%02d_a', pole_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [domain_size_z domain_size_y domain_size_x]);
    h5write (h5_file, dset_name, permute(eval(var_name), [3 2 1]));

    dset_name = sprintf('/materials/debye_%02d_b', pole_index);
    var_name = sprintf('debye_%02d_b', pole_index);
    if exist(var_name, 'var') == 0; 
        error([input_file ': matrix ' var_name ' does not exsit']); 
    end
    h5create(h5_file, dset_name, [domain_size_z domain_size_y domain_size_x]);
    h5write (h5_file, dset_name, permute(eval(var_name), [3 2 1]));
end
end %end of function

