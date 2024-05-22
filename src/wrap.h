#ifndef NEW_WRAPH_H
#define NEW_WRAPH_H

// #include "wrap.h"
#include <array>
#include <exception>
#include <memory>
#include <string>
#include <algorithm>
#include <type_traits>
#include <utility>

#include <TGraphErrors.h>
#include <TH1.h>

namespace ris2{

#ifdef USE_QNTOOLS
using Correlation = Qn::DataContainerStatCalculate;
#endif
using Histogram = TH1;
using Graph = TGraphErrors;

class CannotOpenAFile : public std::exception{ 
public:
  CannotOpenAFile( std::string file_name ){ 
    err_message_.append( file_name );
  }
  CannotOpenAFile( const CannotOpenAFile& ) = default;
  CannotOpenAFile( CannotOpenAFile&& ) = default;
  CannotOpenAFile& operator=( const CannotOpenAFile& ) = default;
  CannotOpenAFile& operator=( CannotOpenAFile&& ) = default;
  ~CannotOpenAFile() = default;
  const char* what() const noexcept override { return err_message_.c_str(); }
private:
  std::string err_message_{"Cannot open a file "};
};

class CannotPullAnObject : public std::exception{ 
public:
  CannotPullAnObject( std::string file_name, std::string object ){ 
    err_message_.append( "Cannot pull an object " ).append( object ).append( " from the file " ).append(file_name);
  }
  CannotPullAnObject( const CannotPullAnObject& ) = default;
  CannotPullAnObject( CannotPullAnObject&& ) = default;
  CannotPullAnObject& operator=( const CannotPullAnObject& ) = default;
  CannotPullAnObject& operator=( CannotPullAnObject&& ) = default;
  ~CannotPullAnObject() = default;
  const char* what() const noexcept override { return err_message_.c_str(); }
private:
  std::string err_message_{};
};

/** 
* @brief Style class is used to set and store information about object draw style
**/
struct Style{
  Style() = default;
  ~Style() = default;
  Style( const Style& ) = default;
  Style( Style&& ) noexcept = default;
  Style& operator=( const Style& ) = default;
  Style& operator=( Style&& ) noexcept = default;
  /// @param color color of the object
  Style& SetColor( int color ){ color_ = color; return *this; }
  /// @param marker is the marker syle of the object; positive is for marker style, negative is for line style 
  Style& SetMarker( int marker ){ marker_ = marker; return *this; }
  Style& operator()( TGraphErrors* points) {
    points->SetLineColor(color_);
    points->SetMarkerColor(color_);
    if( marker_ >= 0 )
      points->SetMarkerStyle(marker_);
    if( marker_ < 0 )
      points->SetLineStyle( abs(marker_) );
    return *this;
  }
  int color_{kBlack};
  int marker_{kFullCircle};
};

/// @brief Dummy class for operating and storing the drawable object and pull the result in the form of TGraph object
/// @param T is a type of drawable object. Supported for TH1*, Qn::DataContainerStatCalculate and TGraph
template<typename T> 
class Result{
public:
  Result() = default;
  template<typename... Args>
  Result( Args... args ){}
  ~Result() = default;
  /// @brief Function for extracting the result points in the form of TGraph. Specialized for each type of storing objects
  TGraphErrors* GetPoints(){ return nullptr; }
  template<typename... Args>
  /// @brief Function for rebinning the underlying object. Specialized for each type of storing objects
  Result<T>& Rebin(Args... args){ return *this; }
  template<typename... Args>
  /// @brief Function for projecting object on one of the axes. Specialized for each type of storing objects
  Result<T>& Project(Args... args){ return *this; }
  /// @brief Function for scaling object. Specialized for each type of storing objects
  Result<T>& Scale( double scale ){ return *this; }
  /// @brief Function used for building ratios. Specialized for each type of storing objects
  Result<T> Divide( const Result<T> other ) const { return Result<T>{}; }
  Result<T> ScaleXaxis( double scale ) const { return Result<T>{}; }
};

/// @brief Dummy class for calculating the systematical variation
/// @param T is a type of drawable object. Supported for TH1*, Qn::DataContainerStatCalculate
template<typename T> 
class Systematics{ 
public:
  Systematics() = default;
  template<typename... Args>
  Systematics( Args... args ){}
  ~Systematics() = default;
  template<typename... Args>
  Systematics<T>& AddToSystematics( Args... args ){ return *this; }
  TGraphErrors* GetPoints(){ return nullptr; }
  template<typename... Args>
  Systematics<T>& Rebin(Args... args){ return *this; }
  template<typename... Args>
  Systematics<T>& Project(Args... args){ return *this; }
  Systematics<T>& Scale( double scale ){ return *this; }
  Systematics<T> ScaleXaxis( double scale ) const { return Systematics<T>{}; }
};

/// @brief Interface class for manipulating the drawable object.
/// @param T is a type of drawable object. Supported for TH1*, Qn::DataContainerStatCalculate and TGraphErrors*
template<typename T>
class Wrap{
public:
  /// @brief A default constructor with no arguments
  Wrap<T>() = default;
  /// @brief A constuctor initializing the underlying object
  /// @param title is the title of the data points
  /// @param Args for Correlation: std::string file_name, std::vector<std::string>> objects, std::vector<double> weights
  /// @param Args for Histogram: std::string file_name, std::vector<std::string>> objects
  /// @param Args for Graph: std::string file_name, std::string> object
  template<typename... Args>
  Wrap<T>( std::string title, Args... args ) : title_{title}, result_{ args... }, systematics_(args...) {  }
  /// @brief Contructor wrapping the pointer to the object. Written for Wrap<Graph> initialization with custom graph.
  Wrap<T>( T* obj ) : result_{ obj }, systematics_(obj) {  }
  /// @brief default destructor
  ~Wrap<T>() = default;
  /// @brief A copy constructor. Copies the underlying objects.
  Wrap<T>( const Wrap<T>& other ) : 
    result_(other.result_), 
    systematics_(other.systematics_),
    style_(other.style_),
    title_(other.title_) {  }
  /// @brief A copy assignement operator. Copies the underlying objects.
  Wrap<T> operator=( const Wrap<T>& other ) { 
    result_ = other.result_;
    systematics_ = other.systematics_;
    style_ = other.style_;
    title_ = other.title_;
    result_points_.reset();
    sys_error_points_.reset();
    return *this;
  }
  /// @brief A move constructor. Moves the underlying objects.
  Wrap<T>( Wrap<T>&& other ) noexcept : 
    result_(std::move(other.result_)), 
    systematics_(std::move(other.systematics_)),
    style_(std::move(other.style_)), 
    title_(std::move(other.title_)) {  }
  /// @brief A move assignement operator. Moves the underlying objects.
  Wrap<T> operator=( Wrap<T>&& other ) noexcept { 
    result_ = std::move(other.result_);
    systematics_ = std::move(other.systematics_);
    style_ = std::move(other.style_);
    title_ = std::move(other.title_);
    result_points_.reset();
    sys_error_points_.reset();
    return *this;
  }
  /// @brief Returns the result in the form TGraphErrors but keeps owning the object.
  TGraphErrors* GetResult() { UpdatePoints(); return result_points_.get(); }
  /// @brief Releases the result points. Memory management is handled to the caller.
  TGraphErrors* ReleaseResult(){ UpdatePoints(); return result_points_.release(); }
  /// @brief Returns the systematics in the form TGraphErrors but keeps owning the object.
  TGraphErrors* GetSystematics() { UpdatePoints(); return sys_error_points_.get(); }
  /// @brief Releases the systematic points. Memory management is handled to the caller.
  TGraphErrors* ReleaseSystematics(){ UpdatePoints(); return sys_error_points_.release(); }
  /// @brief For delayed initialization of the Wrap
  Wrap<T>& SetResult( Result<T>&& res ){ result_ = std::move(res); result_points_.reset(); return *this; };
  /// @brief For delayed initialization of the Wrap
  Wrap<T>& SetSystematics( Systematics<T> sys ){ systematics_ = std::move(sys); sys_error_points_.reset(); };
  const Style& GetStyle(){ return style_; }
  Wrap<T>& SetStyle( Style style ){ style_ = std::move(style); return *this; }
  const std::string& GetTitle() const { return title_; }
  Wrap<T>& SetTitle( std::string title ){ title_ = std::move(title); return *this; }
  template<typename... Args> 
  /// @brief Performs the rebinning of the underlying object. Supported only for Correlation and Histogram
  /// @param Args has std::vector<Qn::AxisD> type for Correlation
  /// @param Args has size_t type for Correlation
  Wrap<T>& Rebin( Args... args ){
    result_.Rebin( args... );
    systematics_.Rebin( args... );
    return *this;
  }
  /// @brief Performs the rebinning of the underlying object. Supported only for Correlation.
  /// @param Args has std::vector<Qn::AxisD> type for Correlation
  template<typename... Args> 
  Wrap<T>& Project( Args... args ){
    result_.Project( args... );
    systematics_.Project( args... );
    return *this;
  }
  /// @brief Scaling the underlying objects. Supported for all three types.
  Wrap<T>& Scale( double scale ){
    result_.Scale(scale);
    systematics_.Scale(scale);
    return *this;
  }
  /// @brief Scaling the X-axis. Supported for all three types.
  Wrap<T>& ScaleXaxis( double scale ){
    result_.ScaleXaxis(scale);
    systematics_.ScaleXaxis(scale);
    return *this;
  }
  /// @brief Fitting the result points.
  Wrap<T>& Fit( TF1* fit_function, std::vector<double> fit_range = {} ){
    UpdatePoints();
    if( fit_range.empty() )
      result_points_->Fit( fit_function, "B", "" );
    else
      result_points_->Fit( fit_function, "B", "", fit_range.at(0), fit_range.at(1) );
    auto res_fit = dynamic_cast<TF1*>( result_points_->GetListOfFunctions()->First() );
    res_fit->SetLineColor(style_.color_);
    return *this;
  }
  Wrap<T>& SetOption( std::string option ){ option_ = std::move(option); return *this; }
  const std::string GetOption(){ return option_; }
private:
  void UpdatePoints(){
    if( !result_points_ ){
      result_points_ = std::unique_ptr<TGraphErrors>( result_.GetPoints() );
    }
    if( !sys_error_points_ ){
      sys_error_points_ = std::unique_ptr<TGraphErrors>( systematics_.GetPoints() );
    }
    if( result_points_ ){
      style_( result_points_.get() );
      result_points_->SetFillColorAlpha(style_.color_, 0.1);
    }
    if( sys_error_points_ ){
      style_( sys_error_points_.get() );
      sys_error_points_->SetFillColorAlpha(style_.color_, 0.1);
    }
  }
  std::string title_{};
  Style style_{};
  Result<T> result_;
  Systematics<T> systematics_;
  std::unique_ptr<TGraphErrors> result_points_{};
  std::unique_ptr<TGraphErrors> sys_error_points_{};
  std::string option_;
};

#ifdef USE_QNTOOLS

template<>
class Result<Qn::DataContainerStatCalculate>{
public:
  Result<Qn::DataContainerStatCalculate>( std::string str_file_name, std::vector<std::string> objects, std::vector<double> weights = {} ){
    auto file_in = std::make_unique<TFile>( str_file_name.c_str(), "READ" );
    if( !file_in )
      throw CannotOpenAFile( str_file_name );
    Qn::DataContainerStatCalculate* ptr_stat_calculate{nullptr};
    Qn::DataContainerStatCollect* ptr_stat_collect{nullptr};

    std::queue<Qn::DataContainerStatCalculate> correlations;
    int i=0;
    for( auto name : objects ){
      auto weight = !weights.empty() ? weights.at(i) : 1.0;
      file_in->GetObject( name.c_str(), ptr_stat_calculate );
      file_in->GetObject( name.c_str(), ptr_stat_collect );
      if( ptr_stat_calculate ){
        correlations.push( Qn::DataContainerStatCalculate(*ptr_stat_calculate)*weight );
        i++;
        continue;
      }
      if( ptr_stat_collect ){
        correlations.push( Qn::DataContainerStatCalculate(*ptr_stat_collect)*weight );
        i++;
        continue;
      }
      throw CannotPullAnObject(str_file_name, name);
    }
    average_ = correlations.front();
    correlations.pop();
    auto list_merge = new TList;
    while( !correlations.empty() ) {
      auto to_merge = new Qn::DataContainerStatCalculate( correlations.front() );
      list_merge->Add(to_merge);
      correlations.pop();
    }
    average_.Merge(list_merge);
  }
  Result<Qn::DataContainerStatCalculate>( const Result<Qn::DataContainerStatCalculate>& ) = default;
  Result<Qn::DataContainerStatCalculate>( Result<Qn::DataContainerStatCalculate>&& ) noexcept = default;
  Result<Qn::DataContainerStatCalculate>& operator=( const Result<Qn::DataContainerStatCalculate>& ) = default;
  Result<Qn::DataContainerStatCalculate>& operator=( Result<Qn::DataContainerStatCalculate>&& ) noexcept = default;
  Result<Qn::DataContainerStatCalculate>( const Qn::DataContainerStatCalculate& container ) : average_{ container } {}

  Result<Qn::DataContainerStatCalculate>& Rebin( std::vector<Qn::AxisD> axes ){
    for( const auto& axis : axes ){
      average_ = average_.Rebin(axis);
    }
    return *this;
  }
  Result<Qn::DataContainerStatCalculate>& Project( std::vector<Qn::AxisD> axes ){
    Rebin( axes );
    auto projection_axes = std::vector<std::string>{};
    for( const auto& axis : axes ){
      projection_axes.emplace_back( axis.Name() );
    }
    average_ = average_.Projection(projection_axes);
    return *this;
  }
  Result<Qn::DataContainerStatCalculate>& Scale(double scale){
    average_ = average_ * scale;
    return *this;
  }
  TGraphErrors* GetPoints() {
    return Qn::ToTGraph( average_ );
  }
  Result<Qn::DataContainerStatCalculate> Divide( Result<Qn::DataContainerStatCalculate> other){
    return Result{ average_ / other.average_ };
  }
private:
  Qn::DataContainerStatCalculate  average_{};
};

template<>
class Systematics<Qn::DataContainerStatCalculate>{
public:
  Systematics<Qn::DataContainerStatCalculate>( std::string str_file_name, std::vector<std::string> objects, std::vector<double> weights = {}  ) {
    auto file_in = std::make_unique<TFile>( str_file_name.c_str(), "READ" );
    if( !file_in )
      throw CannotOpenAFile( str_file_name );
    Qn::DataContainerStatCalculate* ptr_stat_calculate{nullptr};
    Qn::DataContainerStatCollect* ptr_stat_collect{nullptr};
    size_t i=0;
    for( auto name : objects ){
      auto weight = !weights.empty() ? weights.at(i) : 1.0;
      file_in->GetObject( name.c_str(), ptr_stat_calculate );
      file_in->GetObject( name.c_str(), ptr_stat_collect );
      if( ptr_stat_calculate ){
        averaging_objects_.push_back( Qn::DataContainerStatCalculate(*ptr_stat_calculate)*weight );
        ++i;
        continue;
      }
      if( ptr_stat_collect ){
        averaging_objects_.push_back( Qn::DataContainerStatCalculate(*ptr_stat_collect)*weight );
        ++i;
        continue;
      }
      throw CannotPullAnObject(str_file_name, name);
    }
  }
  Systematics<Qn::DataContainerStatCalculate>( const Systematics<Qn::DataContainerStatCalculate>& ) = default;
  Systematics<Qn::DataContainerStatCalculate>( Systematics<Qn::DataContainerStatCalculate>&& ) noexcept = default;
  Systematics<Qn::DataContainerStatCalculate>& operator=( const Systematics<Qn::DataContainerStatCalculate>& ) = default;
  Systematics<Qn::DataContainerStatCalculate>& operator=( Systematics<Qn::DataContainerStatCalculate>&& ) noexcept = default;
  Systematics<Qn::DataContainerStatCalculate>& Rebin( std::vector<Qn::AxisD> axes ){
    for( auto& axis : axes ){
      std::for_each( averaging_objects_.begin(), averaging_objects_.end(), 
        [axis]( auto& obj ){ 
          obj = obj.Rebin(axis); 
        } );
    }
    return *this;
  }
  Systematics<Qn::DataContainerStatCalculate>& Project( std::vector<Qn::AxisD> axes ){
    Rebin( axes );
    auto projection_axes = std::vector<std::string>{};
    for( const auto& axis : axes ){
      projection_axes.emplace_back( axis.Name() );
    }
    std::for_each( averaging_objects_.begin(), averaging_objects_.end(), 
      [projection_axes]( auto& obj ){ 
        obj = obj.Projection(projection_axes); 
      } );
    return *this;
  }
  Systematics<Qn::DataContainerStatCalculate>& Scale(double scale){
    std::for_each( averaging_objects_.begin(), averaging_objects_.end(), 
                  [scale]( auto& obj ){ obj = obj * scale; } );
    return *this;
  }
  TGraphErrors* GetPoints(){
    auto average = averaging_objects_.front();
    auto list_merge = new TList;
    for( int i=1; i<averaging_objects_.size(); ++i ) {
      auto to_merge = new Qn::DataContainerStatCalculate( averaging_objects_.at(i) );
      list_merge->Add(to_merge);
    }
    average.Merge(list_merge);
    std::vector<Qn::DataContainerStatCalculate> variations{};
    std::for_each( averaging_objects_.begin(), averaging_objects_.end(), 
      [&variations, &average] ( Qn::DataContainerStatCalculate& obj ) mutable {
        auto var = average - obj;
        variations.push_back( var );
      } );
    auto systematic_points = Qn::ToTGraph ( average );
    for( int i = 0; i < systematic_points->GetN(); ++i ){
      auto x_hi = average.GetAxes().front().GetUpperBinEdge(i);
      auto x_lo = average.GetAxes().front().GetLowerBinEdge(i);
      auto x_err = fabs( x_hi - x_lo ) / 6;
      auto y_err = fabs(variations.front().At(i).Mean());
      for( const auto& cont : variations ){
        auto err = fabs( cont.At(i).Mean() );
        y_err = std::max( y_err, err );
      }
      systematic_points->SetPointError( i, x_err, y_err );
    }
    return systematic_points;
  }
private:
  std::vector<Qn::DataContainerStatCalculate>  averaging_objects_{};
};

#endif

template<>
class Result<TH1>{
public:
  Result<TH1>( std::string str_file_name, std::vector<std::string> vec_objects ){
    auto file = std::make_unique<TFile>( str_file_name.c_str(), "READ" );
    std::queue<TH1*> histograms;
    TH1* ptr;
    for( const auto& name : vec_objects ){
      file->GetObject( name.c_str(), ptr );
      if( !ptr )
        throw std::runtime_error( std::string("No object in a file").append( " " ).append( name ) );
      histograms.push( ptr );
    }
    auto new_name = std::string(histograms.front()->GetName()) + "_copy";
    histogram_ = std::unique_ptr<TH1>( dynamic_cast<TH1*>(histograms.front()->Clone( new_name.c_str() )) );
    histograms.pop();
    auto list_merge = new TList;
    while( !histograms.empty() ) {
      auto to_merge = histograms.front();
      list_merge->Add(to_merge);
      histograms.pop();
    }
    histogram_->Merge( list_merge );
  }
  Result<TH1>( const Result<TH1>& other ){
    auto new_name = std::string(other.histogram_->GetName()).append("_copy");
    histogram_.reset( dynamic_cast<TH1*>(other.histogram_->Clone( new_name.c_str() )) ); 
  }
  Result<TH1>& operator=( const Result<TH1>& other ) {
    auto new_name = std::string(other.histogram_->GetName()).append("_copy");
    histogram_.reset( dynamic_cast<TH1*>(other.histogram_->Clone( new_name.c_str() )) );
    return *this;
  }
  Result<TH1>( Result<TH1>&& ) = default;
  Result<TH1>& operator=( Result<TH1>&& ) = default;

  TGraphErrors* GetPoints() const {
    auto n_bins = histogram_->GetNbinsX();
    auto x_axis = std::vector<double>{};
    auto y_axis = std::vector<double>{};
    auto y_error = std::vector<double>{};
    for( int i=0; i<n_bins; ++i ){
      x_axis.push_back( histogram_->GetBinCenter(i+1) );
      y_axis.push_back( histogram_->GetBinContent(i+1) );
      y_error.push_back( histogram_->GetBinError(i+1) );
    }
    auto graph = new TGraphErrors(n_bins, x_axis.data(), y_axis.data(), nullptr, y_axis.data() );
    return graph;
  }
  Result<TH1>& Rebin(int n_groups){
    histogram_->RebinX(n_groups);
    return *this;
  }
  Result<TH1>& Scale( double scale ){
    histogram_->Scale(scale);
    return *this;
  }
  Result<TH1> Divide( const Result<TH1>& other ){
    auto result = Result<TH1>(*this);
    result.histogram_->Divide( other.histogram_.get() );
    return result;
  }
private:
  std::unique_ptr<TH1> histogram_;
};

template<>
class Systematics<TH1>{
public:
  Systematics<TH1>( std::string str_file_name, std::vector<std::string> vec_objects ){
    auto file = std::make_unique<TFile>( str_file_name.c_str(), "READ" );
    std::queue<TH1*> histograms;
    TH1* ptr;
    for( const auto& name : vec_objects ){
      file->GetObject( name.c_str(), ptr );
      if( !ptr )
        throw std::runtime_error( std::string("No object in a file").append( " " ).append( name ) );
      averaging_objects_.emplace_back( ptr );
    }
  }
  Systematics<TH1>& Rebin(int n_groups){
    std::for_each( averaging_objects_.begin(), averaging_objects_.end(),
      [n_groups]( const auto& ptr ){ ptr->Rebin(n_groups); } );
    return *this;
  }
  Systematics<TH1>& Scale( double scale ){
    std::for_each( averaging_objects_.begin(), averaging_objects_.end(),
      [scale]( const auto& ptr ){ ptr->Scale(scale); } );
    return *this;
  }
  TGraphErrors* GetPoints(  ){
    auto new_name = std::string(averaging_objects_.front()->GetName()) + "_copy";
    auto average = std::unique_ptr<TH1>( dynamic_cast<TH1*>(averaging_objects_.front()->Clone( new_name.c_str() )) );
    auto list_merge = new TList;
    for( int i=1; i<averaging_objects_.size(); ++i ) {
      auto to_merge = averaging_objects_.at(i).get();
      list_merge->Add(to_merge);
    }
    average->Merge( list_merge );
    auto n_bins = average->GetNbinsX();
    auto x_axis = std::vector<double>{};
    auto y_axis = std::vector<double>{};
    auto y_error = std::vector<double>{};
    for( int i=0; i<n_bins; ++i ){
      x_axis.push_back( average->GetBinCenter(i+1) );
      y_axis.push_back( average->GetBinContent(i+1) );
      y_error.push_back( average->GetBinError(i+1) );
    }
    auto graph = new TGraphErrors(n_bins, x_axis.data(), y_axis.data(), nullptr, y_axis.data() );
    for( int i=0; i<graph->GetN(); ++i ){
      auto x_hi = average->GetXaxis()->GetBinUpEdge(i);
      auto x_lo = average->GetXaxis()->GetBinLowEdge(i);
      auto x_err = fabs(x_hi - x_lo) / 6;
      auto y_err = fabs(average->GetBinContent(i) - averaging_objects_.front()->GetBinContent(i));
      for( const auto& obj : averaging_objects_ ){
        auto err = fabs(average->GetBinContent(i) - obj->GetBinContent(i));
        y_err = std::max( y_err, err );
      }
      graph->SetPointError( i, x_err, y_err );
    }
    return graph;
  }
private:
  std::vector<std::unique_ptr<TH1>> averaging_objects_{};
};

template<>
class Result<TGraphErrors>{
public:
  Result() = default;
  Result( std::string str_file_name, std::string str_obj ){
    auto file = std::make_unique<TFile>( str_file_name.c_str(), "READ" );
    TGraph* ptr{nullptr};
    file->GetObject(str_obj.c_str(), ptr);
    graph_ = std::unique_ptr<TGraphErrors>( dynamic_cast<TGraphErrors*>(ptr->Clone()) );
  }
  Result( TGraphErrors* points ) : graph_( std::unique_ptr<TGraphErrors>(points) ) {  }
  Result( Result<TH1>&  histogram ) : graph_{ std::unique_ptr<TGraphErrors>(  histogram.GetPoints() ) }{ }
  #ifdef QNTOOLS
  Result( Result<Qn::DataContainerStatCalculate>&  correlation ) : graph_{ std::unique_ptr<TGraphErrors>(  correlation.GetPoints() ) }{ }
  #endif
  Result& Scale(double scale){
    for( size_t i=0; i<graph_->GetN(); ++i ){
      auto x = graph_->GetPointX(i);
      auto y = graph_->GetPointY(i); 
      auto y_err = graph_->GetErrorY(i); 
      graph_->SetPoint( i, x,  y*scale );
      graph_->SetPointError( i, 0,  y_err*scale );
    }
    return *this;
  }
  Result& ScaleXaxis(double scale){
    for( size_t i=0; i<graph_->GetN(); ++i ){
      auto x = graph_->GetPointX(i);
      auto y = graph_->GetPointY(i); 
      auto y_err = graph_->GetErrorY(i); 
      graph_->SetPoint( i, x*scale, y );
      graph_->SetPointError( i, 0, y_err );
    }
    return *this;
  }
  Result& SetXAxis( std::vector<double> x_axis ){
    int i=0;
    for( auto x : x_axis ){
      graph_->SetPointX( i, x);
      ++i;
    }
    return *this;
  }
  template<typename Func>
  Result& Perform( const Func& function ){ function(graph_.get()); return *this; }
  TGraphErrors* GetPoints(){
    return dynamic_cast<TGraphErrors*>( graph_->Clone() );
  }
private:
  std::unique_ptr<TGraphErrors> graph_{};
};

}


#endif // NEW_WRAPH_H