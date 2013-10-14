##****************************************************************************
## Time measurements of blocks (imperative style)
## Calling style: tic/toc
## tic (re-)starts time measurement on a specific counter
## toc ends and optionally print time, and further optionally can restart (tic)
## That means that with toc, you can measure time on several timepoints (call
## tic once and toc several times).
## You can use any amount of counters, therefore i.e. time measurements of
## inner and outer loops at the same time is possible. Without any spcification
## a default timer is used.
##
## Examples:
##   tic()
##   ...
##   timepoint[1] <- toc()    ## Time from tic()
##   ...
##   timepoint[2] <- toc()    ## Time from tic()
##
##   tic("loop")              ## Measurement for timer "loop" starting
##   ...
##   tic("another")           ## Measurement for another timer
##   ...
##   toc("loop",echo=TRUE)    ## Print time from tic("loop") to now
##   ...
##   toc("another",echo=TRUE) ## Print time from toc("another") till now
##
##   toc("loop", echo=TRUE, text="Time is passing..." )
##       ## This prints: "Time is passing... 2d 07:52:42"  (with your time of course)
##****************************************************************************

.time.temp <- list( default=NULL )

## Start time measurement (or reset it)
tic <- function( timer="default" ) {
  .time.temp[[timer]] <<- proc.time()[3]
}

## Tick-Tock: return time span from last tic. Optionally print it.
## Parameters: timer     - Name of the timer.
##             tic       - Call tick at the end of toc
##             echo      - Print time
##             echo.full - Human readable time format
##             text      - Optional text printed before the time (if echo=TRUE)
toc <- function( timer="default", tic=FALSE, echo=FALSE, echo.full=TRUE, text="" ) {
  ## If tic() was never called.
  if( is.null( .time.temp[[timer]] ) ) {
    warning( paste( "tic not called before toc. Timer:", timer ) )
    tic(timer)
  }

  time <- proc.time()[3] - .time.temp[[timer]]

  if( echo ) {
    ## Full format is hh:mm::ss or D hh::mm::ss if running for more than one day.
    if( echo.full ) {
      hours <- as.integer( floor( time / (60*60)) %% 60 )

      print( sprintf( '%s%s%s:%02d:%02d',
                      text,
                      ifelse( text != "", " ", "" )[1],
                      ifelse( hours >= 24,
                              sprintf( "%dd %02d",
                                      as.integer(floor(hours/24)),
                                      as.integer(hours %% 24)), hours ),
                      as.integer( floor( time / 60 ) %% 60 ),
                      as.integer( time %% 60 )
            ) )
    } else {
      print( paste( text, time ) )  ## Simple time print.
    }

    if( tic )
      tic(timer)

    return( invisible( time ) )     ## Printed, so invisible return
  }
  else {
    if( tic )
      tic(timer)

    return( time )
  }
}

## Alias functions for tic() and toc()
tick <- function( ... ) tic( ... )
tack <- function( ... ) toc( ... )

## Just print the current time in the specific timer -- without modifying it.
## Uses toc()
timer <- function( timer="default" ) {
  toc( timer=timer, tic=FALSE, echo=TRUE, echo.full=TRUE )
}
