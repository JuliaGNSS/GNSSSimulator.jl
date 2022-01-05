struct Order{N} end
Order(N) = Order{N}()

function get_process(::Order{1}, T)
    1
end

function get_process(::Order{2}, T)
    @SMatrix [
        1 T
        0 1
    ]
end

function get_process(::Order{3}, T)
    @SMatrix [
        1 T T^2 / 2
        0 1 T
        0 0 1
    ]
end

function get_process_covariance(::Order{1}, T)
    @SMatrix [T] # SMatrix more efficient than scalar for cholesky
end

function get_process_covariance(::Order{2}, T)
    @SMatrix [
        T^3/3 T^2/2
        T^2/2 T
    ]
end

function get_process_covariance(::Order{3}, T)
    @SMatrix [
        T^5/20 T^4/8 T^3/6
        T^4/8 T^3/3 T^2/2
        T^3/6 T^2/2 T
    ]
end