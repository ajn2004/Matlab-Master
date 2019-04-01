function [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc(i2, thrds, rads, lpcnt)

        [m,~] = size(i2);
        switch m % Switch localization algorithm based on size

            case 9
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_3(i2, thrds, rads, lpcnt);
            case 25
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_5(i2, thrds, rads, lpcnt);
            case 49
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_7(i2, thrds, rads, lpcnt);
            case 81
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_9(i2, thrds, rads, lpcnt);
            case 121
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_11(i2, thrds, rads, lpcnt);
            case 169
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_13(i2, thrds, rads, lpcnt);
            case 225
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_15(i2, thrds, rads, lpcnt);   
            case 529
                [xf, xc, yf, yc, Np, Nc, sx, sxc, sy, syc, off, offc, lv] = slim_chain_loc_23(i2, thrds, rads, lpcnt);   
            otherwise
                xf = [];
                xc = [];
                yf = [];
                yc = [];
                Np = [];
                Nc = [];
                sx = [];
                sxc = [];
                sy = [];
                syc = [];
                off = [];
                offc = [];
                lv = [];
        end