function g = GravityForce(mass)
% Gravity force: @ -mg

    g = (-1).*mass.*PhysicsConstants.GravityConst;
end