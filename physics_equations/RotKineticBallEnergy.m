function r = RotKineticBallEnergy(mass, circ_speed, radius)
% Ball rotation energy: @ E = Iw^2
    r = (2/5).*mass.*(radius^2).*(circ_speed.^2);
end
