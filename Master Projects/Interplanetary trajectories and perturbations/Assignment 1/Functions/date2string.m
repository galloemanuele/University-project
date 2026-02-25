function str = date2string(dateVec)
% DESCRIPTION:
% Function to convert the specified date to a string
% ------------------------------------------------------------------
% PROTOTYPE:
% str = date2string(dateVec)
% -------------------------------------------------------------------
% INPUT:
%  dateVec[6x1]     Date in the Gregorian calendar
%                   [year, month, day, hour, minute, second]
% -------------------------------------------------------------------
% OUTPUT:
%  str[1]           String with the name of the month corresponding to the
%                   second component of the date vector.
% -------------------------------------------------------------------
% CONTRIBUTORS:
%  Alessia Casiero
%  Elisa Fiume 
%  Emanuele Gallo 
%  Arianna Marotta
% -------------------------------------------------------------------
% VERSIONS:
%  01-01-2024: First version
% -------------------------------------------------------------------
months = {'January', 'February', 'March', 'April', 'May', 'June', ...
	'July', 'August', 'September', 'October', 'November', 'December'};
str = sprintf('%i %s %i at %02i:%02i:%02i', dateVec(1), months{dateVec(2)}, dateVec(3), dateVec(4), dateVec(5), round(dateVec(6)));

end

