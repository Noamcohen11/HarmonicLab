function h = HeightEnergy(mass, height)
% Height energy: @ E = mgh
    h = mass.*height.*PhysicsConstants.GravityConst;
end