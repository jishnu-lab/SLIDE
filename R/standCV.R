#' Standardization for k-fold cross-validation.
#'
#' Standardize the training and validation sets for cross-validation.
#'
#' @param valid_y the validation set made from the response vector
#' @param valid_x the validation set made from the data matrix
#' @param train_y the training set made from the response vector
#' @param train_x the training set made from the data matrix
#' @export

standCV <- function(valid_x, train_x, valid_y = NULL, train_y = NULL)  {
    train_x_stand <- scale(train_x, T, T)
    centers_x <- attr(train_x_stand, "scaled:center")
    scales_x <- attr(train_x_stand, "scaled:scale")
    valid_x_stand <- t((t(valid_x) - centers_x) / scales_x)

    if (!is.null(valid_y)) {
        train_y_stand <- scale(train_y, T, T)
        centers_y <- attr(train_y_stand, "scaled:center")
        scales_y <- attr(train_y_stand, "scaled:scale")
        valid_y_stand <- t((t(valid_y) - centers_y) / scales_y)
        return (list("valid_y" = valid_y_stand,
                     "valid_x" = valid_x_stand,
                     "train_y" = train_y_stand,
                     "train_x" = train_x_stand))
    }

    return (list("valid_x" = valid_x_stand,
                 "train_x" = train_x_stand))
}
