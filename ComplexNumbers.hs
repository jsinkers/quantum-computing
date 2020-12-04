module ComplexNumbers
  (
    ComplexNumber,
    Re,
    Im,
    Rad,
    Theta,
    ComplexVector,
    ComplexMatrix,
    complexAdd,
    complexMult,
    complexSub,
    complexModulus,
    complexConjugate,
    complexInverse,
    complexPower,
    complexRoots,
    complexDiv,
    complexRe,
    complexIm,
    complexRadius,
    complexPhase,
    toPolar,
    toCartesian,
    cvAdd,
    cvScalarMult,
    cvInverse,
    cvDot,
    cmIdentity,
    cmAdd,
    cmScalarMult,
    cmInverse,
    cmTranspose,
    cmDagger,
    cmAdjoint,
    cmConjugate,
    cmMult,
    printComplexMatrix,
    showComplexMatrix
   ) where 

import Data.List (transpose, intercalate)
-- data types
-- ============================================================================
type Re = Double
type Im = Double

type Rad = Double
type Theta = Double

data ComplexNumber = Cartesian Re Im | Euler Rad Theta

type ComplexVector = [ComplexNumber]

-- m x n complex matrix: contains m lists of n elements
type ComplexMatrix = [[ComplexNumber]]

instance Show ComplexNumber where
    show (Cartesian r i) = show r ++ " + " ++ show i ++ "i"
    show (Euler r t) = show r ++ "e^{i" ++ show t ++ "}"

-- complex number algebra
-- ============================================================================
-- complex zero number
complexZero :: ComplexNumber
complexZero = Cartesian 0 0

complexAdd :: ComplexNumber -> ComplexNumber -> ComplexNumber
complexAdd c1 c2 = Cartesian (r1+r2) (i1+i2)
    where (Cartesian r1 i1) = toCartesian c1 
          (Cartesian r2 i2) = toCartesian c2

complexMult :: ComplexNumber -> ComplexNumber -> ComplexNumber
complexMult (Euler r1 t1) (Euler r2 t2) = Euler (r1*r2) (t1+t2)
complexMult c1 c2 = Cartesian (r1*r2 - i1*i2) (r1*i2 + r2*i1)
    where (Cartesian r1 i1) = toCartesian c1
          (Cartesian r2 i2) = toCartesian c2

complexSub :: ComplexNumber -> ComplexNumber -> ComplexNumber
complexSub c1 c2 = Cartesian (r1-r2) (i1-i2)
    where (Cartesian r1 i1) = toCartesian c1
          (Cartesian r2 i2) = toCartesian c2

complexModulus :: ComplexNumber -> Double
complexModulus (Cartesian r i) = sqrt (r^2 + i^2)
complexModulus (Euler r _) = r

-- determine complex conjugate of a complex number
complexConjugate :: ComplexNumber -> ComplexNumber
complexConjugate (Cartesian r1 i1) = Cartesian r1 (-i1)
complexConjugate (Euler r t) = Euler r (-t)

-- determine inverse of a complex number (i.e. z^-1 = 1/z) - not to be confused with additive inverse
complexInverse :: ComplexNumber -> ComplexNumber
complexInverse c = complexMult (Cartesian (1/modCsq) 0) (complexConjugate c)
    where modCsq = complexModulus c ^ 2

-- division of complex numbers
complexDiv :: ComplexNumber -> ComplexNumber -> ComplexNumber
complexDiv c1 c2 = complexMult c1 $ complexInverse c2

-- convert to polar form
toPolar :: ComplexNumber -> ComplexNumber
toPolar p@(Euler _ _) = p
toPolar c@(Cartesian r i) = Euler (complexModulus c) (atan2 i r)

-- convert to Cartesian form
toCartesian :: ComplexNumber -> ComplexNumber
toCartesian c@(Cartesian _ _) = c
toCartesian (Euler r theta) = Cartesian (r*cos theta) (r*sin theta)

-- getter for real part of a complex number
complexRe :: ComplexNumber -> Re
complexRe c = r
    where (Cartesian r _) = toCartesian c

-- getter for imaginary part of a complex number
complexIm :: ComplexNumber -> Im
complexIm c = i
    where (Cartesian _ i) = toCartesian c

-- getter for magnitude of a complex number
complexRadius :: ComplexNumber -> Rad
complexRadius = complexModulus

-- getter for phase of a complex number
complexPhase :: ComplexNumber -> Theta
complexPhase c = t
    where (Euler _ t) = toPolar c

complexPower :: Double -> ComplexNumber -> ComplexNumber
complexPower n c = Euler (r**n) (n*t)
    where (Euler r t) = toPolar c

-- need a function to convert angles to (-Pi, Pi] or [0,2Pi]


-- determine roots of unity
complexRoots :: Int -> [ComplexNumber]
complexRoots n = [Euler 1 (fromIntegral (2*k) * pi/fromIntegral n) | k <- [0..n-1]]

-- complex vector space operations
-- ============================================================================
cvAdd :: ComplexVector -> ComplexVector -> ComplexVector
cvAdd cv1s cv2s = if length cv1s == length cv2s 
                     then [complexAdd c1 c2 | (c1, c2) <- zip cv1s cv2s]
                     else error "ComplexVector size mismatch"

cvScalarMult :: ComplexNumber -> ComplexVector -> ComplexVector
cvScalarMult s = map (complexMult s) 

cvInverse :: ComplexVector -> ComplexVector
cvInverse = map (complexSub (Euler 0 0)) 

-- dot product of two complex vectors
cvDot :: ComplexVector -> ComplexVector -> ComplexNumber
cvDot cv1 cv2 =  foldl complexAdd complexZero $ zipWith complexMult cv1 cv2

-- complex matrix functions
-- ============================================================================
-- complex identity matrix I_n
cmIdentity :: Int -> ComplexMatrix
cmIdentity n = [[if i == j then Cartesian 1 0 else complexZero | j <- [0..n-1]] | i <- [0..n-1]]

-- determine result of addition of 2 complex matrices
cmAdd :: ComplexMatrix -> ComplexMatrix -> ComplexMatrix
cmAdd cm1s cm2s = if length cm1s == length cm2s 
                     then [cvAdd c1 c2 | (c1, c2) <- zip cm1s cm2s]
                     else error "ComplexMatrix size mismatch"

-- determine scalar multiple of a complex matrix
cmScalarMult :: ComplexNumber -> ComplexMatrix -> ComplexMatrix
cmScalarMult s = map (cvScalarMult s)

-- determine additive inverse of a complex matrix
cmInverse :: ComplexMatrix -> ComplexMatrix
cmInverse = map cvInverse

-- determine transpose of a complex matrix
cmTranspose :: ComplexMatrix -> ComplexMatrix
cmTranspose = transpose

-- determine conjugate of a complex matrix
cmConjugate :: ComplexMatrix -> ComplexMatrix
cmConjugate = map (map complexConjugate) 

-- determine adjoint/dagger
cmDagger :: ComplexMatrix -> ComplexMatrix
cmDagger = cmConjugate . cmTranspose

cmAdjoint :: ComplexMatrix -> ComplexMatrix
cmAdjoint = cmDagger

showComplexMatrix :: ComplexMatrix -> String
showComplexMatrix cm = intercalate "\n" $ map show cm 

printComplexMatrix :: ComplexMatrix -> IO ()
printComplexMatrix cm = do putStr $ showComplexMatrix cm
                           putStrLn ""

-- multiply 2 complex matrices
-- note no size checking is performed
cmMult :: ComplexMatrix -> ComplexMatrix -> ComplexMatrix
cmMult cm1 cm2 = [[cvDot cv1 cv2 | cv2 <- cm2'] | cv1 <- cm1]
    where cm2' = cmTranspose cm2

-- Action of application of complex matrix A to complex vector v
cmAction :: ComplexMatrix -> ComplexVector -> ComplexVector
cmAction a v = [cvDot ra v | ra <- a]

-- test matrix
cm1 :: ComplexMatrix
cm1 = [[Cartesian 6 (-3), Cartesian 2 12,  Cartesian 0 (-19) ],
                 [Cartesian 0 0,    Cartesian 5 2.1, Cartesian 17 0    ],
                 [Cartesian 1 0,    Cartesian 2 5,   Cartesian 3 (-4.5)]]

cm2 :: ComplexMatrix
cm2 = [[Cartesian 1 0, Cartesian 0 0,  Cartesian 0 0 ],
                 [Cartesian 0 0,    Cartesian 1 0.0, Cartesian 0 0    ],
                 [Cartesian 0 0,    Cartesian 0 0,   Cartesian 1 0]]

cm3 :: ComplexMatrix
cm3 = [[Cartesian 1 0, Cartesian 1 1],
       [Cartesian 1 1, Cartesian 1 0]]

cm4 :: ComplexMatrix
cm4 = [[Cartesian 2 1, Cartesian 4 1],
       [Cartesian 3 0, Cartesian 5 0]]

