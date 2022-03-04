function v = convo_PC(u1,u2,alpha_mat)

if not( size(u1) == size(u2))
    error("u1 and u2 should have the same size")
end

v = vect2tens( convo_PC_mat(u1,alpha_mat)*tens2vect(u2), size(u1) );

    