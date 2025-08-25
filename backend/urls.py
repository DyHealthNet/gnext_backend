from .views.variant_view import PheWASView, VariantAnnotationView
from .views.overview_view import OverviewDataView
from .views.trait_view import ManhattanView, QQView, TraitView, ChromosomeBoundsView
from django.urls import path


urlpatterns = [
    path("variant_phewas/", PheWASView.as_view(), name="variant_phewas"),
    path("variant_annotation/", VariantAnnotationView.as_view(), name="variant_annotation"),
    path("overview_stats/", OverviewDataView.as_view(), name="overview_stats"),
    path("trait_manhattan/", ManhattanView.as_view(), name="trait_manhattan"),
    path("trait_qq/", QQView.as_view(), name="trait_qq"),
    path("trait_get_variants/", TraitView.as_view(), name="trait_get_variants"),
    path("trait_get_chromosomeBounds/", ChromosomeBoundsView.as_view(), name="trait_get_chromosomeBounds"),
]

