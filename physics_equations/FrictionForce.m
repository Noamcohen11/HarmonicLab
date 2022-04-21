function f = FrictionForce(friction_const, normal)
%% Friction force = -const.*normal.
    f = friction_const.*normal;
end