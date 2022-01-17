function filtered_signal = notch_filt(signal, fs, f0)
    filtered_signal = signal;
    n = 10; %filter order
    bw = 3; %bandwidth
    for i = 1 : length(f0)
        notchSpecs = fdesign.notch('N,F0,BW', n, f0(i), bw, fs);
        notchFilt = design(notchSpecs, 'butter', 'Systemobject', true);
        filtered_signal = notchFilt(filtered_signal);
    end
end

