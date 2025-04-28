from django.urls import path
from api.views.typesense_autocomplete import AutoComplete

urlpatterns = [
    path('autocomplete/', AutoComplete.as_view()),
]