{-# OPTIONS_GHC -Wall #-}

{-# LANGUAGE NegativeLiterals #-}

import           Numeric.Sundials.CVode.ODE
import           Numeric.Sundials.ODEOpts (ODEOpts(..))
import           Numeric.LinearAlgebra


roberts :: Double -> Vector Double -> Vector Double
roberts t v = vector $ robertsAux t (toList v)
  where
    robertsAux _ [y1, y2, y3] = [ -0.04 * y1 + 1.0e4 * y2 * y3
                         , 0.04 * y1 - 1.0e4 * y2 * y3 - 3.0e7 * (y2)^(2 :: Int)
                         , 3.0e7 * (y2)^(2 :: Int)
                         ]
    robertsAux _ _ = error "roberts RHS not defined"

rootFn :: Double -> Vector Double -> Vector Double
rootFn _ v = case xs of
               [y1, _y2, y3] -> vector [ y1 - 0.0001
                                       , y3 - 0.01
                                       ]
               _             -> error "roberts root function RHS not defined"
  where
    xs = toList v

ts :: [Double]
ts = take 12 $ map (* 10.0) (0.4 : ts)

solve :: SolverResult Matrix Vector Int Double
solve = odeSolveRootVWith' opts BDF (ScXX' 1.0 1.0e-4 1.0 1.0 (vector [1.0e-8, 1.0e-14, 1.0e-6]))
                      Nothing roberts (vector [1.0, 0.0, 0.0])
                      2 rootFn
                      (vector ts)
  where
    opts = ODEOpts { maxNumSteps = 10000
                   , minStep     = 1.0e-12
                   , relTol      = 1.0e-4
                   , absTols     = vector [1.0e-8, 1.0e-14, 1.0e-6]
                   , initStep    = Nothing
                   , maxFail     = 10
                   }

main :: IO ()
main = putStrLn $ show solve

