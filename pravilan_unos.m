function y = pravilan_unos(x,xmin,xmax)

switch nargin
    
    case 1
        if (isempty(x))
            fprintf('\n GRESKA: Vrednost nije uneta! \n');
            y = 0;
        elseif (~isnumeric(x))
            fprintf('\n GRESKA: Potrebno je uneti numericku vrednost! \n');
            y = 0;
        else
            y = 1;
        end
        
    case 2
        if (isempty(x))
            fprintf('\n GRESKA: Vrednost nije uneta! \n');
            y = 0;
        elseif (~isnumeric(x))
            fprintf('\n GRESKA: Potrebno je uneti numericku vrednost! \n');
            y = 0;
        elseif (x < xmin)
            fprintf('\n GRESKA: Uneta vrednost nije u zadatom opsegu! \n');
            fprintf('\n Uneta vrednost mora biti veca od ili jednaka  %.2f. \n', xmin);
            y = 0;
        else
            y = 1;
        end
        
    case 3
     
        if (isempty(x))
            fprintf('\n GRESKA: Vrednost nije uneta! \n');
            y = 0;
        elseif (~isnumeric(x))
            fprintf('\n GRESKA: Potrebno je uneti numericku vrednost! \n');
            y = 0;
        elseif (x < xmin || x > xmax)
            fprintf('\n GRESKA: Uneta vrednost nije u zadatom opsegu! \n');
            fprintf('\n Uneta vrednost mora biti veca od ili jednaka %.2f. \n', xmin);
            fprintf('\n Uneta vrednost mora biti manja od ili jednaka %.2f. \n', xmax);
            y = 0;
        else
            y = 1;
        end
        
end

