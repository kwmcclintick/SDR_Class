% AUTHOR: Kyle, Nivetha, Krysta
% Date: 2/20/18
% Description: BBMAC_TX takes in a frame that is ready for transmission
% and attaches a CRC to it, as well as the index of the best carrier in
% a list of carriers to use for transmission based on energy states.
% Carriers are assumed to be ordered by ascending frequency.

function [best_carrier, energy_states, frame_with_crc] = BBMAC_TX(energy_levels, energy_states, FC, SRC, DEST, PAYLOAD, t)
% ADD CRC BITS TO FRAME
gen = comm.CRCGenerator([1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1],'ChecksumsPerFrame',1);
frame_to_send = FC;
frame_to_send = [frame_to_send SRC];
frame_to_send = [frame_to_send DEST];
frame_to_send = [frame_to_send PAYLOAD];
frame_with_crc = step(gen,frame_to_send');

% MEMORYLESS BUMBLEBEE MODEL
R = 10e3; C = 100e-6;
Tau = R*C; % energy state decay mimics a capacitor's voltage over 1 second

if energy_levels > energy_states
    energy_states = energy_levels;
else
    energy_states = energy_states + ...
        (energy_states - energy_levels)*(1 - t/exp(1/Tau));
end
[min_e_carrier, min_e_carrier_idx] = min(energy_states);
best_carrier = min_e_carrier_idx;
end

