function pix2pho = em_gain(emg)

switch emg
    case 0
        pix2pho = 1;
    case 50
        pix2pho = 6.94;
    case 100
        pix2pho = 12.87;
    case 150
        pix2pho = 18.42;
    case 200
        pix2pho = 23.56;
    case 250
        pix2pho = 28.42;
    case 300
        pix2pho = 33.33;
    otherwise
        pix2pho = [];
end