fixes_transition_fit <- function() {
    paths <- as.matrix(expand.grid(rep(list(1:2), 4L)))
    people <- rep(seq_len(nrow(paths)), each = 10L)
    residual <- rep(stats::qnorm((seq_len(10L) - 0.5)/10), nrow(paths))
    data <- data.frame(id = rep(seq_along(people), each = 4L), 
        time = rep(seq_len(4L), length(people)), y = 0.1 * (as.vector(t(2 * 
            paths[people, ] - 3)) + 0.5 * rep(residual, each = 4L)))
    lta(data, "y", "id", n_profiles = 2, time = "time", transitions = "occasion", 
        n_starts = 1, seed = 1)
}
test_that("transition parameter_inference() uses the caller's step in its table", 
    {
        skip_on_cran()
        fit <- fixes_transition_fit()
        expect_s3_class(fit, "multilpa_lta")
        expect_true(fit$converged)
        expect_false(fit$boundary)
        expect_false(fit$probability_boundary)
        expect_lt(max(abs(stats::plogis(fit$initial_coefficients) - 
            0.5)), 0.25)
        expect_lt(max(abs(stats::plogis(fit$transition_coefficients) - 
            0.5)), 0.25)
        inference <- .lta_inference
        steps <- numeric()
        local_mocked_bindings(.lta_inference = function(fit, 
            vcov_type, step = 1e-04) {
            steps <<- c(steps, step)
            inference(fit, vcov_type, step)
        })
        coarse <- parameter_inference(fit, step = 0.01)
        fine <- parameter_inference(fit, step = 1e-04)
        expect_identical(steps, c(0.01, 1e-04))
        covariance <- attr(coarse, "covariance_unconstrained")
        expect_equal(coarse$standard_error[coarse$block == "transition"], 
            unname(sqrt(diag(covariance)))[startsWith(colnames(covariance), 
                "transition.")], tolerance = 1e-12)
        expect_equal(coarse$standard_error, fine$standard_error, 
            tolerance = 1e-07)
    })
