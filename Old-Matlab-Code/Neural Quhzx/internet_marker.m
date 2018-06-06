email_address = {'andrew.j.nelson@maine.edu','desiree.runka@maine.edu'};
[reg, name] = system('hostname');
message = ['Your internet is up'];
subject = ['INTERNET FOR YOU'];
while true
    try
        matlabmail(email_address, message, subject , 'matlab4fpalm@gmail.com','Fpalm4life');
    catch lsterr
    end
end