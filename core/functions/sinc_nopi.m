function f = sinc_nopi(x)
    % Modified sinc function because matlab sinc function already includes the factor pi.
    % This makes notation consistent with definitions in the publications.
    f = sinc(x/pi);
end

