function s = import2d(name, pop_time, time_stamp, undersampling, zeropad_times, mess_date, scans)

s = load2d(name, pop_time, time_stamp, 'meta', true, 'noise', true, 'scans', scans);
s.date = mess_date;
s.undersampling = undersampling;
s.zeropad = zeropad_times * length(s.time);
s = freq2d(s);
s.w3 = s.freq - 0;
s.phase = s.phase + pi;