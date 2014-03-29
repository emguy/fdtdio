function status = fdtd_run(input_file, output_file)
% NOTE: This program is free software; you can redistribute it
% and/or modify it under the terms of the GNU General Public 
% License as published by the Free Software Foundation; either
% version 3, or (at you option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Pulic License for more details.

full_path =  mfilename('fullpath');
dir = full_path(1:end - length(mfilename()));
if strcmp(computer(), 'PCWIN64'); 
    excutable = [dir 'ezfdtd.exe'];
else
    excutable = [dir 'ezfdtd'];
end
lock_file = [input_file(1:end-2) 'lock'];
status = 0;

display(excutable);
% ---------------------------------------
% Check files
% ---------------------------------------
if exist(output_file, 'file') == 2;
    delete(output_file); % remove the file to avoid overwrite
end
if exist(lock_file, 'file') == 2;
    delete(lock_file); % remove the file to avoid overwrite
end
fclose(fopen(output_file, 'w')); % create an empty file
lock = fopen(lock_file, 'w'); fclose(lock);

% ---------------------------------------
% simulate on a separate console
% ---------------------------------------
% windows platform 
if strcmp(computer(), 'PCWIN64'); 
    command_before = ['start cmd /k "title "' input_file, '" && '];
    command = [excutable ' ' input_file ' ' output_file];
    command_after = [' && del ' lock_file];
    command_close = ['&& exit "'];
    system([command_before command command_after command_close]);

    % wait for the simulation
    while exist(lock_file)==2; %#ok<EXIST>
        pause(2);
    end

else if strcmp(computer(), 'GLNXA64'); 
    command = [excutable ' ' input_file ' ' output_file ' && rm -f ' lock_file];
    system(command);
else
    disp('Incompatible platmorm. Exit !');
    return;
end

status = 1;

end
