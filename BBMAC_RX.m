% AUTHOR: Kyle, Nivetha, Krysta
% Date: 2/20/18
% Description: BBMAC_RX confirms theres no CRC error in a frame given to it
% by the physical layer team, makes sure the radio receiving the frame has
% the same SRC ID as the DEST ID of the incoming frame, and divides the
% frame up into subsections for processing by other teams

function [FC, SRC, DEST, PAYLOAD, CRC, overhead_flag] = BBMAC_RX(frame, our_SRC)
overhead_flag = false;

% check to see if there is a crc error
detect = comm.CRCDetector([1 1 0 0 0 0 0 0 0 0 0 0 0 0 1 1],'ChecksumsPerFrame',1);
[message, CRC] = step(detect,frame);
frame = message;
if CRC == 0
    DEST = frame(13:20);
    if DEST == our_SRC
        % Frame divided according to task 3's proposal
        FC = frame(1:4);
        SRC = frame(5:12);
        PAYLOAD = frame(21:end);
    else
        % if destination is wrong, DEST returned as all zeros
        DEST = zeros(1, 8);
    end
end

% check for overhead message
if SRC == DEST
    overhead_flag = true;
end
end

