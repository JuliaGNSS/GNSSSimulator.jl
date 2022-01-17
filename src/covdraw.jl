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

function soft_bound(value, noise, upper_bound, lower_bound, relative_soft_bounding = 0.2)
    value += value[1] > upper_bound * (1 - sign(upper_bound) * relative_soft_bounding) && value[end] > 0 && noise[end] > 0 ?
        SVector(noise[1], -noise[end]) : noise
    value += value[1] < lower_bound * (1 + sign(lower_bound) * relative_soft_bounding) && value[end] < 0 && noise[end] < 0 ?
        SVector(noise[1], -noise[end]) : noise
    value
end