function i2 = sub_image(i1,cents,pixw)

wind = -pixw:pixw;
i2 = i1(cents(2) + wind, cents(1) + wind);
