using System;
using System.Collections.Generic;
using System.Text;

namespace UMF_lab2;

public interface IBasisFunction
{
    double Value(double x);
    double Derivative(double x);
}

public interface IFiniteElement
{
    int Id { get; }
    IReadOnlyList<int> GlobalNodes { get; }
    IReadOnlyList<IBasisFunction> BasisFunctions { get; }
    double Length { get; }
}

// Базовый интерфейс для всех сеток
public interface IMesh
{
    // Количество узлов/точек в сетке
    int NodeCount { get; }
}

// Пространственная сетка
public interface ISpaceMesh : IMesh
{
    IReadOnlyList<double> Nodes { get; }
    IReadOnlyList<IFiniteElement> Elements { get; }
    int ElementCount { get; }
    
    // Левая граница области
    double LeftBoundary { get; }
    
    // Правая граница области
    double RightBoundary { get; }

    double GetNodeCoordinate(int nodeId);
    IFiniteElement GetElement(int elementId);
}

// Временная сетка
public interface ITimeMesh : IMesh
{
    // Все временные точки: t₀, t₁, ..., tₙ
    IReadOnlyList<double> TimePoints { get; }
    
    // Количество временных шагов (точек - 1)
    int TimeStepCount { get; }
    
    // Начальное время
    double TimeStart { get; }
    
    // Конечное время
    double TimeEnd { get; }
    
    // Получить шаг по времени между слоями n и n+1
    double GetTimeStep(int stepIndex);

    // Получить время в точке n
    double GetTimePoint(int stepIndex);
}

public interface ICoefficients
{
    double Lambda(double u, double x, double t = 0);
    double Sigma(double u, double x, double dudx, double t = 0);
    double SigmaDuDt(double u, double x, double t = 0);
    double F(double u, double x, double t = 0);
    double Theta(double u, double x, double t = 0);
}

public interface IBoundaryConditions
{
    BoundaryType LeftType { get; }
    BoundaryType RightType { get; }
    double LeftValue { get; }
    double RightValue { get; }
    Func<double, double, double> Theta { get; }
    bool IsBoundaryNode(int nodeId, int totalNodes);
}

public enum BoundaryType
{
    Dirichlet,   // 1-го рода
    Neumann,     // 2-го рода
    Robin        // 3-го рода
}

public interface ILocalAssembly
{
    double[,] LocalMatrix { get; }
    double[] LocalVector { get; }
    int Size { get; }

    void Compute(IFiniteElement element, ICoefficients coeffs,
                    IReadOnlyList<double> globalSolution);
}

public interface IGlobalAssembly
{
    BandedMatrix GlobalMatrix { get; }
    double[] GlobalVector { get; }

    void Assemble(ISpaceMesh mesh, ILocalAssembly localAssembly,
        ICoefficients coeffs, IReadOnlyList<double> solution,
        IBoundaryConditions bc);

    void ApplyBoundaryConditions(IBoundaryConditions bc, IReadOnlyList<double> solution);
}

public interface ILinearSolver
{
    double[] Solve(BandedMatrix matrix, double[] rhs);
}

public interface INonlinearSolver
{
    NonlinearSolverResult Solve(INonlinearProblem problem, double[] initialGuess);
}

public interface INonlinearProblem
{
    int Size { get; }
    void ComputeSystem(IReadOnlyList<double> q, out BandedMatrix A, out double[] b);
    double[] ComputeResidual(IReadOnlyList<double> q);
    double ComputeResidualNorm(IReadOnlyList<double> q);
}

public interface ISolutionResult
{
    double[] Solution { get; }
    double[] ResidualHistory { get; }
    int IterationsCount { get; }
    bool Converged { get; }
    double FinalResidual { get; }
    double ComputationTime { get; }
}

public interface IInitialCondition
{
    // Начальное условие u(0, x) = g_0(x)
    double Value(double x);
}

public interface ITimeIntegrator
{
    // Выполняет один временной шаг
    double[] TimeStep(INonlinearProblem problem, double[] uCurrent, double dt, double time);

    // Возвращает историю невязок за последний шаг
    double[] LastResidualHistory { get; }
    
    // Количество итераций за последний временной шаг
    int LastIterationsCount { get; }
}

public interface IParabolicProblem : INonlinearProblem
{
    // Коэффициент при du/dt в параболическом уравнении
    double GetMassCoefficient(double u, double x, double t);
    
    // Вычисляет систему для параболического уравнения
    void ComputeParabolicSystem(IReadOnlyList<double> uPrev, double dt, double time,
        out BandedMatrix A, out double[] b);
}

public interface IElementFactory
{
    IFiniteElement CreateElement(int id, int node1, int node2);
    IBasisFunction CreateLinearBasisFunction(int localId, double x1, double x2);
    IBasisFunction CreateQuadraticBasisFunction(int localId, double x1, double x2, double x3);
}
